from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
import psycopg2
import psycopg2.extras
import os
import urllib.parse
from typing import Optional
import json
import math
import threading

app = FastAPI(title="Wildebeest Migration & SOC Analysis API — Dissertation Edition")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Railway injects DATABASE_URL — parse it if present, else fall back to individual vars
DATABASE_URL = os.getenv("DATABASE_URL")
if DATABASE_URL:
    r = urllib.parse.urlparse(DATABASE_URL)
    DB_CONFIG = {
        "host": r.hostname,
        "port": r.port or 5432,
        "database": r.path.lstrip("/"),
        "user": r.username,
        "password": r.password,
    }
else:
    DB_CONFIG = {
        "host": os.getenv("DB_HOST", "localhost"),
        "port": os.getenv("DB_PORT", "5432"),
        "database": os.getenv("DB_NAME", "MaraSerengeti"),
        "user": os.getenv("DB_USER", "postgres"),
        "password": os.getenv("DB_PASSWORD", "1149"),
    }

@app.get("/")
def root():
    return FileResponse("index.html")

# Migration geom is stored as UTM Zone 36S but labelled SRID=4326 — correct it:
MIGRATION_TO_WGS84 = "ST_Transform(ST_SetSRID(geom, 32736), 4326)"

SEASONS = {
    "janmar":  {"label": "Jan–Mar", "table_suffix": "janmar2019",  "months": [1, 2, 3]},
    "aprmay":  {"label": "Apr–May", "table_suffix": "aprmay2019",  "months": [4, 5]},
    "junjuly": {"label": "Jun–Jul", "table_suffix": "junjuly2019", "months": [6, 7]},
    "augoct":  {"label": "Aug–Oct", "table_suffix": "augoct2019",  "months": [8, 9, 10]},
    "novdec":  {"label": "Nov–Dec", "table_suffix": "novdec2019",  "months": [11, 12]},
}

NDVI_TABLES = {
    "janmar":  "ndvi_janmar_2019",
    "aprmay":  "ndvi_aprmay_2019",
    "junjuly": "ndvi_junjuly_2019",
    "augoct":  "ndvi_augoct_2019",
    "novdec":  "ndvi_novdec_2019",
}

MONTHLY_RAINFALL = {
    1: 61.3, 2: 57.0, 3: 89.2, 4: 186.7, 5: 64.3,
    6: 49.0, 7: 10.8, 8: 36.0, 9: 25.8, 10: 49.2,
    11: 133.8, 12: 233.2
}

SEASON_RAINFALL = {
    "janmar":  round(sum(MONTHLY_RAINFALL[m] for m in [1,2,3]) / 3, 1),
    "aprmay":  round(sum(MONTHLY_RAINFALL[m] for m in [4,5]) / 2, 1),
    "junjuly": round(sum(MONTHLY_RAINFALL[m] for m in [6,7]) / 2, 1),
    "augoct":  round(sum(MONTHLY_RAINFALL[m] for m in [8,9,10]) / 3, 1),
    "novdec":  round(sum(MONTHLY_RAINFALL[m] for m in [11,12]) / 2, 1),
}

NEXT_SEASON = {
    "janmar": "aprmay", "aprmay": "junjuly", "junjuly": "augoct",
    "augoct": "novdec", "novdec": "janmar",
}

# Bounding box for the Mara-Serengeti study area
BBOX = "lon_bin BETWEEN 29 AND 42 AND lat_bin BETWEEN -7 AND 5"
BBOX_XY = "ST_X(geom) BETWEEN 29 AND 42 AND ST_Y(geom) BETWEEN -7 AND 5"


def get_conn(timeout_ms=25000):
    """Connect with statement_timeout so slow queries fail fast instead of hanging."""
    try:
        return psycopg2.connect(**DB_CONFIG, options=f"-c statement_timeout={timeout_ms}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"DB connection failed: {str(e)}")

def dict_cursor(conn):
    return conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

def pearson(xs, ys):
    n = len(xs)
    if n < 2: return None
    mx, my = sum(xs)/n, sum(ys)/n
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    dx = math.sqrt(sum((x-mx)**2 for x in xs))
    dy = math.sqrt(sum((y-my)**2 for y in ys))
    if dx == 0 or dy == 0: return None
    return round(num/(dx*dy), 4)

def linear_regression(xs, ys):
    n = len(xs)
    if n < 2: return None, None, None
    mx, my = sum(xs)/n, sum(ys)/n
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    denom = sum((x-mx)**2 for x in xs)
    if denom == 0: return None, None, None
    slope = num/denom
    intercept = my - slope*mx
    ss_res = sum((y-(slope*x+intercept))**2 for x,y in zip(xs,ys))
    ss_tot = sum((y-my)**2 for y in ys)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
    return round(slope,4), round(intercept,4), round(r2,4)


# ── STARTUP: preload KDE janmar + protected areas only ────────────────────────
# Everything else loads on demand when the user toggles that layer.
# Raster queries use TABLESAMPLE (not ORDER BY random) — much faster.
_cache = {"kde_janmar": None, "protected_areas": None}

def _raw_conn():
    return psycopg2.connect(**DB_CONFIG, options="-c statement_timeout=25000")

def _warm_startup():
    def _run():
        print("[startup] Preloading KDE janmar + protected areas…")

        # KDE janmar — pure vector aggregation, no rasters, fast
        try:
            months = SEASONS["janmar"]["months"]
            mf = f"AND EXTRACT(MONTH FROM TO_DATE(date, 'DD/MM/YYYY')) = ANY(ARRAY{months}::int[])"
            conn = _raw_conn()
            cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
            cur.execute(f"""
                WITH pts AS (
                    SELECT ST_X({MIGRATION_TO_WGS84}) AS lon,
                           ST_Y({MIGRATION_TO_WGS84}) AS lat
                    FROM migration_data WHERE geom IS NOT NULL {mf}
                ),
                binned AS (
                    SELECT ROUND(lon::numeric/0.05)*0.05 AS lon_bin,
                           ROUND(lat::numeric/0.05)*0.05 AS lat_bin,
                           COUNT(*) AS cnt
                    FROM pts WHERE lon BETWEEN 29 AND 42 AND lat BETWEEN -7 AND 5
                    GROUP BY lon_bin, lat_bin
                )
                SELECT lon_bin AS lon, lat_bin AS lat, cnt
                FROM binned ORDER BY cnt DESC LIMIT 8000
            """)
            rows = [dict(r) for r in cur.fetchall()]
            conn.close()
            nonzero = sorted(int(r["cnt"]) for r in rows if int(r["cnt"]) > 0)
            n = len(nonzero)
            p33 = nonzero[int(n*0.33)] if n >= 3 else (nonzero[0] if nonzero else 1)
            p66 = nonzero[int(n*0.66)] if n >= 3 else (nonzero[-1] if nonzero else 2)
            result = []
            for r in rows:
                c = int(r["cnt"])
                cls = 1 if c==0 else (2 if c<=p33 else (3 if c<=p66 else 4))
                result.append({"lat": float(r["lat"]), "lon": float(r["lon"]), "count": c, "class": cls})
            _cache["kde_janmar"] = {"data": result, "breaks": [0, p33, p66], "season": "janmar"}
            print(f"[startup] KDE janmar: {len(result)} cells")
        except Exception as e:
            print(f"[startup] KDE janmar failed: {e}")

        # Protected areas — try UTM→WGS84 transform, fall back to raw geom
        try:
            conn = _raw_conn()
            cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
            for sql in [
                'SELECT name, protection, status, ST_AsGeoJSON(ST_Transform(ST_SetSRID(geom, 32736), 4326)) AS geojson FROM "mara serengeti protected areas"',
                'SELECT name, protection, status, ST_AsGeoJSON(geom) AS geojson FROM "mara serengeti protected areas"',
            ]:
                try:
                    cur.execute(sql)
                    features = []
                    for row in cur.fetchall():
                        g = json.loads(row["geojson"])
                        # Sanity check: coordinates must be in WGS84 range
                        coords = g.get("coordinates", [])
                        if coords:
                            def first_coord(c):
                                while isinstance(c[0], list): c = c[0]
                                return c
                            c0 = first_coord(coords)
                            if not (-180 <= c0[0] <= 180 and -90 <= c0[1] <= 90):
                                continue
                        features.append({
                            "type": "Feature",
                            "properties": {"name": row["name"], "protection": row["protection"], "status": row["status"]},
                            "geometry": g
                        })
                    if features:
                        _cache["protected_areas"] = {"type": "FeatureCollection", "features": features}
                        print(f"[startup] Protected areas: {len(features)} features")
                        break
                except Exception as ie:
                    print(f"[startup] PA attempt failed: {ie}")
                    continue
            conn.close()
        except Exception as e:
            print(f"[startup] Protected areas failed: {e}")

        print("[startup] Done.")

    threading.Thread(target=_run, daemon=True).start()


@app.on_event("startup")
def startup_event():
    _warm_startup()


# ── HEALTH ─────────────────────────────────────────────────────────────────────
@app.get("/api/health")
def health():
    conn = get_conn(timeout_ms=5000)
    conn.close()
    return {
        "status": "ok",
        "cache": {
            "kde_janmar": len(_cache["kde_janmar"]["data"]) if _cache["kde_janmar"] else 0,
            "protected_areas": len(_cache["protected_areas"]["features"]) if _cache["protected_areas"] else 0,
        }
    }


# ── STUDY AREA ─────────────────────────────────────────────────────────────────
@app.get("/api/study-area")
def study_area():
    conn = get_conn()
    cur = conn.cursor()
    for sql in [
        'SELECT ST_AsGeoJSON(ST_Transform(ST_Union(ST_SetSRID(geom, 32736)), 4326)) FROM "mara serengeti protected areas"',
        'SELECT ST_AsGeoJSON(ST_Transform(ST_Union(geom), 4326)) FROM "mara serengeti protected areas"',
        'SELECT ST_AsGeoJSON(ST_Union(geom)) FROM "mara serengeti protected areas"',
    ]:
        try:
            cur.execute(sql)
            row = cur.fetchone()
            if row and row[0]:
                conn.close()
                return {"geojson": json.loads(row[0])}
        except:
            continue
    conn.close()
    return {"geojson": None}


# ── PROTECTED AREAS ────────────────────────────────────────────────────────────
@app.get("/api/protected-areas")
def protected_areas():
    # Serve from startup cache — instant
    if _cache["protected_areas"]:
        return _cache["protected_areas"]
    # Cache not ready yet — query live
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        for sql in [
            'SELECT name, protection, status, ST_AsGeoJSON(ST_Transform(ST_SetSRID(geom, 32736), 4326)) AS geojson FROM "mara serengeti protected areas"',
            'SELECT name, protection, status, ST_AsGeoJSON(geom) AS geojson FROM "mara serengeti protected areas"',
        ]:
            try:
                cur.execute(sql)
                features = [
                    {"type": "Feature",
                     "properties": {"name": r["name"], "protection": r["protection"], "status": r["status"]},
                     "geometry": json.loads(r["geojson"])}
                    for r in cur.fetchall()
                ]
                if features:
                    return {"type": "FeatureCollection", "features": features}
            except:
                continue
        return {"type": "FeatureCollection", "features": []}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ── MIGRATION KDE CLASSES ──────────────────────────────────────────────────────
@app.get("/api/migration/kde-classes")
def migration_kde_classes(season: Optional[str] = "janmar"):
    if season not in SEASONS:
        raise HTTPException(status_code=400, detail="Invalid season")
    # Janmar served from startup cache — instant
    if season == "janmar" and _cache["kde_janmar"]:
        return _cache["kde_janmar"]
    months = SEASONS[season]["months"]
    mf = f"AND EXTRACT(MONTH FROM TO_DATE(date, 'DD/MM/YYYY')) = ANY(ARRAY{months}::int[])"
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute(f"""
            WITH pts AS (
                SELECT ST_X({MIGRATION_TO_WGS84}) AS lon, ST_Y({MIGRATION_TO_WGS84}) AS lat
                FROM migration_data WHERE geom IS NOT NULL {mf}
            ),
            binned AS (
                SELECT ROUND(lon::numeric/0.05)*0.05 AS lon_bin,
                       ROUND(lat::numeric/0.05)*0.05 AS lat_bin,
                       COUNT(*) AS cnt
                FROM pts WHERE lon BETWEEN 29 AND 42 AND lat BETWEEN -7 AND 5
                GROUP BY lon_bin, lat_bin
            )
            SELECT lon_bin AS lon, lat_bin AS lat, cnt FROM binned ORDER BY cnt DESC LIMIT 8000
        """)
        rows = [dict(r) for r in cur.fetchall()]
        conn.close()
        if not rows: return {"data": [], "breaks": []}
        nonzero = sorted(int(r["cnt"]) for r in rows if int(r["cnt"]) > 0)
        n = len(nonzero)
        p33 = nonzero[int(n*0.33)] if n >= 3 else (nonzero[0] if nonzero else 1)
        p66 = nonzero[int(n*0.66)] if n >= 3 else (nonzero[-1] if nonzero else 2)
        result = []
        for r in rows:
            c = int(r["cnt"])
            cls = 1 if c==0 else (2 if c<=p33 else (3 if c<=p66 else 4))
            result.append({"lat": float(r["lat"]), "lon": float(r["lon"]), "count": c, "class": cls})
        return {"data": result, "breaks": [0, p33, p66], "season": season}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        try: conn.close()
        except: pass


# ── SOC LAYER (on demand) ─────────────────────────────────────────────────────
@app.get("/api/soc/sample")
def soc_sample():
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute("""
            SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS soc_value
            FROM (
                SELECT (ST_PixelAsCentroids(rast)).*
                FROM o_4_soc_dataset TABLESAMPLE SYSTEM(8)
            ) AS pixels
            WHERE val > 0 AND val < 500
            LIMIT 4000
        """)
        return {"data": [dict(r) for r in cur.fetchall()]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ── SOC RESIDUALS (on demand) ─────────────────────────────────────────────────
@app.get("/api/soc-residuals/sample")
def soc_residuals_sample():
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute("""
            SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS residual
            FROM (
                SELECT (ST_PixelAsCentroids(rast)).*
                FROM "o_16_soc_residuals_(idwfrom5000points)" TABLESAMPLE SYSTEM(12)
            ) AS pixels
            WHERE val IS NOT NULL
            LIMIT 6000
        """)
        return {"data": [dict(r) for r in cur.fetchall()]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ── MIGRATION POINTS (on demand) ──────────────────────────────────────────────
@app.get("/api/migration/points")
def migration_points(season: Optional[str] = "janmar"):
    if season not in SEASONS:
        raise HTTPException(status_code=400, detail="Invalid season")
    months = SEASONS[season]["months"]
    mf = f"AND EXTRACT(MONTH FROM TO_DATE(date, 'DD/MM/YYYY')) = ANY(ARRAY{months}::int[])"
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute(f"""
            WITH pts AS (
                SELECT ST_X({MIGRATION_TO_WGS84}) AS lon, ST_Y({MIGRATION_TO_WGS84}) AS lat
                FROM migration_data WHERE geom IS NOT NULL {mf}
            )
            SELECT lon, lat FROM pts
            WHERE lon BETWEEN 29 AND 42 AND lat BETWEEN -7 AND 5
            LIMIT 20000
        """)
        rows = cur.fetchall()
        return {"data": [{"lon": float(r["lon"]), "lat": float(r["lat"])} for r in rows], "count": len(rows)}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ── RAINFALL LAYER (on demand) ────────────────────────────────────────────────
@app.get("/api/rainfall/sample")
def rainfall_sample(season: Optional[str] = "janmar"):
    if season not in SEASONS:
        raise HTTPException(status_code=400, detail="Invalid season")
    table = f"{SEASONS[season]['table_suffix']}Rainfall"
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute(f"""
            SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS rainfall_mm
            FROM (
                SELECT (ST_PixelAsCentroids(rast)).*
                FROM "{table}" TABLESAMPLE SYSTEM(15)
            ) AS pixels
            WHERE val > 0
            LIMIT 3000
        """)
        return {"data": [dict(r) for r in cur.fetchall()]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ── NDVI LAYER (on demand) ────────────────────────────────────────────────────
@app.get("/api/ndvi/sample")
def ndvi_sample(season: Optional[str] = "janmar"):
    table = NDVI_TABLES.get(season, "ndvi_janmar_2019")
    conn = get_conn()
    cur = dict_cursor(conn)
    for tbl in [table, "ndvi_janmar_2019"]:
        try:
            cur.execute(f"""
                SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS ndvi
                FROM (
                    SELECT (ST_PixelAsCentroids(rast)).*
                    FROM {tbl} TABLESAMPLE SYSTEM(10)
                ) AS pixels
                WHERE val BETWEEN -1 AND 1
                LIMIT 3000
            """)
            rows = cur.fetchall()
            conn.close()
            return {"data": [dict(r) for r in rows], "table": tbl}
        except:
            continue
    conn.close()
    return {"data": [], "table": None}


# ── DEM LAYER (on demand) ─────────────────────────────────────────────────────
@app.get("/api/dem/sample")
def dem_sample():
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute("""
            SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS elevation
            FROM (
                SELECT (ST_PixelAsCentroids(rast)).*
                FROM (SELECT rast FROM "30mdem" TABLESAMPLE SYSTEM(5) LIMIT 6) tiles
            ) AS pixels
            WHERE val > 0
            LIMIT 2000
        """)
        return {"data": [dict(r) for r in cur.fetchall()]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ── SEASONAL ANALYSIS ──────────────────────────────────────────────────────────
@app.get("/api/analysis/seasonal")
def seasonal_analysis():
    conn = get_conn()
    cur = dict_cursor(conn)
    results = []
    for season_key, season_info in SEASONS.items():
        months = season_info["months"]
        mf = f"EXTRACT(MONTH FROM TO_DATE(date, 'DD/MM/YYYY')) = ANY(ARRAY{months}::int[])"
        try:
            cur.execute(f"SELECT COUNT(*) AS n FROM migration_data WHERE geom IS NOT NULL AND {mf}")
            n = int(cur.fetchone()["n"])
        except: n = 0
        avg_soc = avg_ndvi = avg_elev = avg_ndvi_lag = None
        try:
            cur.execute(f"""
                WITH sample AS (
                    SELECT {MIGRATION_TO_WGS84} AS gw FROM migration_data
                    WHERE geom IS NOT NULL AND {mf} LIMIT 300
                )
                SELECT AVG((SELECT ST_Value(s.rast, p.gw) FROM o_4_soc_dataset s WHERE ST_Intersects(s.rast, p.gw) LIMIT 1)) AS avg_soc,
                       AVG((SELECT ST_Value(d.rast, p.gw) FROM "30mdem" d WHERE ST_Intersects(d.rast, p.gw) LIMIT 1)) AS avg_elev
                FROM sample p
            """)
            row = cur.fetchone()
            if row:
                avg_soc = float(row["avg_soc"]) if row["avg_soc"] else None
                avg_elev = float(row["avg_elev"]) if row["avg_elev"] else None
        except Exception as e: print(f"SOC/elev ({season_key}): {e}")
        ndvi_tbl = NDVI_TABLES.get(season_key, "ndvi_janmar_2019")
        try:
            cur.execute(f"""
                WITH sample AS (SELECT {MIGRATION_TO_WGS84} AS gw FROM migration_data WHERE geom IS NOT NULL AND {mf} LIMIT 200)
                SELECT AVG((SELECT ST_Value(n.rast, p.gw) FROM {ndvi_tbl} n WHERE ST_Intersects(n.rast, p.gw) LIMIT 1)) AS avg_ndvi FROM sample p
            """)
            r2 = cur.fetchone()
            avg_ndvi = float(r2["avg_ndvi"]) if r2 and r2["avg_ndvi"] else None
        except Exception as e: print(f"NDVI ({season_key}): {e}")
        next_s = NEXT_SEASON[season_key]
        lag_tbl = NDVI_TABLES.get(next_s, "ndvi_janmar_2019")
        try:
            cur.execute(f"""
                WITH sample AS (SELECT {MIGRATION_TO_WGS84} AS gw FROM migration_data WHERE geom IS NOT NULL AND {mf} LIMIT 200)
                SELECT AVG((SELECT ST_Value(n.rast, p.gw) FROM {lag_tbl} n WHERE ST_Intersects(n.rast, p.gw) LIMIT 1)) AS avg_ndvi_lag FROM sample p
            """)
            r3 = cur.fetchone()
            avg_ndvi_lag = float(r3["avg_ndvi_lag"]) if r3 and r3["avg_ndvi_lag"] else None
        except Exception as e: print(f"NDVI lag ({season_key}): {e}")
        results.append({
            "season": season_key, "label": season_info["label"], "n_points": n,
            "avg_soc": round(avg_soc, 2) if avg_soc else None,
            "avg_rainfall": SEASON_RAINFALL[season_key],
            "avg_elevation": round(avg_elev, 1) if avg_elev else None,
            "avg_ndvi": round(avg_ndvi, 4) if avg_ndvi else None,
            "avg_ndvi_lag": round(avg_ndvi_lag, 4) if avg_ndvi_lag else None,
            "next_season_label": SEASONS[next_s]["label"],
        })
    conn.close()
    return {"data": results}


# ── GRAZED vs UNGRAZED SOC ─────────────────────────────────────────────────────
@app.get("/api/analysis/grazed-vs-ungrazed")
def grazed_vs_ungrazed():
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute(f"""
            WITH all_pts AS (
                SELECT ROUND(ST_X({MIGRATION_TO_WGS84})::numeric/0.1)*0.1 AS lon_bin,
                       ROUND(ST_Y({MIGRATION_TO_WGS84})::numeric/0.1)*0.1 AS lat_bin,
                       COUNT(*) AS cnt
                FROM migration_data WHERE geom IS NOT NULL GROUP BY lon_bin, lat_bin
            ),
            grid AS (
                SELECT lon_bin, lat_bin, cnt, ST_SetSRID(ST_MakePoint(lon_bin, lat_bin), 4326) AS pt
                FROM all_pts WHERE lon_bin BETWEEN 29 AND 42 AND lat_bin BETWEEN -7 AND 5
            )
            SELECT g.lon_bin, g.lat_bin, g.cnt,
                (SELECT ST_Value(s.rast, g.pt) FROM o_4_soc_dataset s WHERE ST_Intersects(s.rast, g.pt) LIMIT 1) AS soc,
                (SELECT ST_Value(d.rast, g.pt) FROM "30mdem" d WHERE ST_Intersects(d.rast, g.pt) LIMIT 1) AS elevation
            FROM grid g ORDER BY g.cnt DESC LIMIT 2000
        """)
        rows = [dict(r) for r in cur.fetchall()]
        conn.close()
    except Exception as e:
        conn.close(); raise HTTPException(status_code=500, detail=str(e))

    if not rows: return {"data": [], "summary": {}}
    nonzero = sorted(int(r["cnt"]) for r in rows if int(r["cnt"]) > 0)
    n = len(nonzero)
    p33 = nonzero[int(n*0.33)] if n >= 3 else 1
    p66 = nonzero[int(n*0.66)] if n >= 3 else 2
    valid = []
    for r in rows:
        c = int(r["cnt"])
        soc = float(r["soc"]) if r["soc"] and 0 < float(r["soc"]) < 500 else None
        elev = float(r["elevation"]) if r["elevation"] else None
        if soc is None or elev is None: continue
        cls = 1 if c==0 else (2 if c<=p33 else (3 if c<=p66 else 4))
        valid.append({"lon": float(r["lon_bin"]), "lat": float(r["lat_bin"]), "count": c,
                      "class": cls, "grazed": cls > 1, "soc": soc, "elevation": elev})
    if not valid: return {"data": [], "summary": {}}
    slope, intercept, r2 = linear_regression([v["elevation"] for v in valid], [v["soc"] for v in valid])
    for v in valid:
        v["residual_soc"] = round(v["soc"] - (slope*v["elevation"]+intercept), 4) if slope is not None else 0.0
    def gs(pts):
        if not pts: return {}
        socs=[p["soc"] for p in pts]; resids=[p["residual_soc"] for p in pts]; elevs=[p["elevation"] for p in pts]; n=len(pts)
        return {"n":n, "mean_soc":round(sum(socs)/n,2), "mean_residual_soc":round(sum(resids)/n,2),
                "mean_elevation":round(sum(elevs)/n,1),
                "sd_soc":round(math.sqrt(sum((x-sum(socs)/n)**2 for x in socs)/n),2) if n>1 else 0,
                "sd_residual_soc":round(math.sqrt(sum((x-sum(resids)/n)**2 for x in resids)/n),2) if n>1 else 0}
    lm = {1:"Ungrazed",2:"Mild",3:"Moderate",4:"Heavy"}
    return {"data": valid, "summary": {
        "grazed": gs([v for v in valid if v["grazed"]]),
        "ungrazed": gs([v for v in valid if not v["grazed"]]),
        "by_class": {lm[c]: gs([v for v in valid if v["class"]==c]) for c in [1,2,3,4]},
        "ols_slope": slope, "ols_intercept": intercept, "ols_r2": r2,
        "corr_count_residual_soc": pearson([v["count"] for v in valid],[v["residual_soc"] for v in valid]),
        "corr_count_raw_soc": pearson([v["count"] for v in valid],[v["soc"] for v in valid]),
        "breaks": [0, p33, p66]}}


# ── NDVI GRAZED vs UNGRAZED ────────────────────────────────────────────────────
@app.get("/api/analysis/ndvi-grazing")
def ndvi_grazing():
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute(f"""
            WITH all_pts AS (
                SELECT ROUND(ST_X({MIGRATION_TO_WGS84})::numeric/0.1)*0.1 AS lon_bin,
                       ROUND(ST_Y({MIGRATION_TO_WGS84})::numeric/0.1)*0.1 AS lat_bin,
                       COUNT(*) AS cnt
                FROM migration_data WHERE geom IS NOT NULL GROUP BY lon_bin, lat_bin
            )
            SELECT lon_bin, lat_bin, cnt, ST_SetSRID(ST_MakePoint(lon_bin, lat_bin), 4326) AS pt
            FROM all_pts WHERE lon_bin BETWEEN 29 AND 42 AND lat_bin BETWEEN -7 AND 5
        """)
        cells = [dict(r) for r in cur.fetchall()]
    except Exception as e:
        conn.close(); raise HTTPException(status_code=500, detail=str(e))
    nonzero = sorted(int(c["cnt"]) for c in cells if int(c["cnt"]) > 0)
    n = len(nonzero); p33 = nonzero[int(n*0.33)] if n >= 3 else 1
    grazed_pts = [c for c in cells if int(c["cnt"]) > p33]
    ungrazed_pts = [c for c in cells if int(c["cnt"]) <= p33]
    results = []
    for season_key in SEASONS:
        ndvi_tbl = NDVI_TABLES.get(season_key, "ndvi_janmar_2019")
        lag_tbl = NDVI_TABLES.get(NEXT_SEASON[season_key], "ndvi_janmar_2019")
        def get_ndvi_for_pts(pts, tbl, limit=150):
            if not pts: return None
            sample = pts[:limit]
            pt_list = ",".join(f"ST_SetSRID(ST_MakePoint({p['lon_bin']},{p['lat_bin']}),4326)" for p in sample)
            try:
                cur.execute(f"""
                    WITH pts AS (SELECT unnest(ARRAY[{pt_list}]) AS gw)
                    SELECT AVG((SELECT ST_Value(n.rast, p.gw) FROM {tbl} n WHERE ST_Intersects(n.rast, p.gw) LIMIT 1)) AS avg_ndvi FROM pts p
                """)
                row = cur.fetchone(); return float(row["avg_ndvi"]) if row and row["avg_ndvi"] else None
            except: return None
        gn=get_ndvi_for_pts(grazed_pts,ndvi_tbl); un=get_ndvi_for_pts(ungrazed_pts,ndvi_tbl)
        gl=get_ndvi_for_pts(grazed_pts,lag_tbl); ul=get_ndvi_for_pts(ungrazed_pts,lag_tbl)
        results.append({"season":season_key,"label":SEASONS[season_key]["label"],
            "next_label":SEASONS[NEXT_SEASON[season_key]]["label"],"rainfall":SEASON_RAINFALL[season_key],
            "grazed_ndvi":round(gn,4) if gn else None,"ungrazed_ndvi":round(un,4) if un else None,
            "grazed_ndvi_lag":round(gl,4) if gl else None,"ungrazed_ndvi_lag":round(ul,4) if ul else None})
    conn.close(); return {"data": results}


# ── CORRELATION DATA ───────────────────────────────────────────────────────────
@app.get("/api/analysis/correlation")
def correlation_data():
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute(f"""
            WITH grid AS (
                SELECT ROUND(ST_X({MIGRATION_TO_WGS84})::numeric,1) AS lon_bin,
                       ROUND(ST_Y({MIGRATION_TO_WGS84})::numeric,1) AS lat_bin,
                       COUNT(*) AS grazing_count,
                       ST_SetSRID(ST_MakePoint(AVG(ST_X({MIGRATION_TO_WGS84})),AVG(ST_Y({MIGRATION_TO_WGS84}))),4326) AS center_pt
                FROM migration_data WHERE geom IS NOT NULL GROUP BY lon_bin, lat_bin
            )
            SELECT lon_bin, lat_bin, grazing_count,
                (SELECT ST_Value(s.rast, g.center_pt) FROM o_4_soc_dataset s WHERE ST_Intersects(s.rast, g.center_pt) LIMIT 1) AS soc,
                (SELECT ST_Value(r.rast, g.center_pt) FROM janmar2019Rainfall r WHERE ST_Intersects(r.rast, g.center_pt) LIMIT 1) AS rainfall,
                (SELECT ST_Value(d.rast, g.center_pt) FROM "30mdem" d WHERE ST_Intersects(d.rast, g.center_pt) LIMIT 1) AS elevation,
                (SELECT ST_Value(n.rast, g.center_pt) FROM ndvi_janmar_2019 n WHERE ST_Intersects(n.rast, g.center_pt) LIMIT 1) AS ndvi
            FROM grid g WHERE lon_bin BETWEEN 29 AND 42 AND lat_bin BETWEEN -7 AND 5
            ORDER BY grazing_count DESC LIMIT 800
        """)
        return {"data": [dict(r) for r in cur.fetchall()]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ── DEBUG ──────────────────────────────────────────────────────────────────────
@app.get("/api/debug-env")
def debug_env():
    return {"DATABASE_URL": "SET" if os.environ.get("DATABASE_URL") else "NOT SET"}

@app.get("/api/debug/protected-areas")
def debug_protected_areas():
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute("SELECT column_name FROM information_schema.columns WHERE table_name = 'mara serengeti protected areas' ORDER BY ordinal_position")
        cols = [r['column_name'] for r in cur.fetchall()]
        non_geom = [c for c in cols if c != 'geom']
        col_list = ", ".join(f'"{c}"' for c in non_geom)
        cur.execute(f'SELECT {col_list}, ST_SRID(geom) as srid, ST_XMin(ST_Envelope(geom)) as xmin, ST_YMin(ST_Envelope(geom)) as ymin FROM "mara serengeti protected areas" LIMIT 3')
        return {"columns": cols, "sample": [dict(r) for r in cur.fetchall()], "cache_loaded": bool(_cache["protected_areas"])}
    except Exception as e:
        return {"error": str(e)}
    finally:
        conn.close()

@app.get("/api/debug/tables")
def debug_tables():
    conn = get_conn()
    cur = dict_cursor(conn)
    try:
        cur.execute("SELECT table_name FROM information_schema.tables WHERE table_schema = 'public' ORDER BY table_name")
        return {"all_tables": [r['table_name'] for r in cur.fetchall()], "cache": {k: bool(v) for k,v in _cache.items()}}
    except Exception as e:
        return {"error": str(e)}
    finally:
        conn.close()
        
