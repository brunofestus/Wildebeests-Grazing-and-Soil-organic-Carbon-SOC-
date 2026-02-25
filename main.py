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

app = FastAPI(title="Wildebeest Migration & SOC Analysis API — Dissertation Edition")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

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

# Migration geom is SRID=4326 but contains UTM Zone 36S values — fix:
MIGRATION_TO_WGS84 = "ST_Transform(ST_SetSRID(geom, 32736), 4326)"

SEASONS = {
    "janmar":  {"label": "Jan–Mar", "table_suffix": "janmar2019",  "months": [1, 2, 3]},
    "aprmay":  {"label": "Apr–May", "table_suffix": "aprmay2019",  "months": [4, 5]},
    "junjuly": {"label": "Jun–Jul", "table_suffix": "junjuly2019", "months": [6, 7]},
    "augoct":  {"label": "Aug–Oct", "table_suffix": "augoct2019",  "months": [8, 9, 10]},
    "novdec":  {"label": "Nov–Dec", "table_suffix": "novdec2019",  "months": [11, 12]},
}

# NDVI table for each season (for lagged analysis, we look at T+1 NDVI after T grazing)
NDVI_TABLES = {
    "janmar":  "ndvi_janmar_2019",
    "aprmay":  "ndvi_aprmay_2019",
    "junjuly": "ndvi_junjuly_2019",
    "augoct":  "ndvi_augoct_2019",
    "novdec":  "ndvi_novdec_2019",
}

# Monthly rainfall from dissertation Table 2
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

# Season labels for lagged (T+1) analysis
NEXT_SEASON = {
    "janmar": "aprmay",
    "aprmay": "junjuly",
    "junjuly": "augoct",
    "augoct": "novdec",
    "novdec": "janmar",
}


def get_conn():
    try:
        return psycopg2.connect(**DB_CONFIG)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"DB connection failed: {str(e)}")

def dict_cursor(conn):
    return conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return None
    mx = sum(xs)/n
    my = sum(ys)/n
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    dx = math.sqrt(sum((x-mx)**2 for x in xs))
    dy = math.sqrt(sum((y-my)**2 for y in ys))
    if dx == 0 or dy == 0:
        return None
    return round(num/(dx*dy), 4)

def linear_regression(xs, ys):
    """Returns slope, intercept, r-squared"""
    n = len(xs)
    if n < 2:
        return None, None, None
    mx = sum(xs)/n
    my = sum(ys)/n
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    denom = sum((x-mx)**2 for x in xs)
    if denom == 0:
        return None, None, None
    slope = num/denom
    intercept = my - slope*mx
    ss_res = sum((y - (slope*x+intercept))**2 for x,y in zip(xs,ys))
    ss_tot = sum((y-my)**2 for y in ys)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
    return round(slope,4), round(intercept,4), round(r2,4)

# ─── HEALTH ────────────────────────────────────────────────────────────────────
@app.get("/api/health")
def health():
    conn = get_conn()
    conn.close()
    return {"status": "ok"}


# ─── STUDY AREA ────────────────────────────────────────────────────────────────
@app.get("/api/study-area")
def study_area():
    conn = get_conn()
    cur = conn.cursor()
    for sql in [
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


# ─── MIGRATION DENSITY — 4 GRAZING CLASSES (Jenks-like via percentile breaks) ─
@app.get("/api/migration/kde-classes")
def migration_kde_classes(season: Optional[str] = "janmar"):
    """
    Returns binned GPS points per 0.05° cell, then classifies into:
    Class 1 = ungrazed (0 points), Class 2 = mild, 3 = moderate, 4 = heavy
    using percentile breaks on non-zero counts (matching dissertation KDE + Jenks approach).
    """
    if season not in SEASONS:
        raise HTTPException(status_code=400, detail="Invalid season")
    months = SEASONS[season]["months"]
    month_filter = f"AND EXTRACT(MONTH FROM TO_DATE(date, 'DD/MM/YYYY')) = ANY(ARRAY{months}::int[])"

    conn = get_conn()
    cur = dict_cursor(conn)
    query = f"""
        WITH pts AS (
            SELECT
                ST_X({MIGRATION_TO_WGS84}) AS lon,
                ST_Y({MIGRATION_TO_WGS84}) AS lat
            FROM migration_data
            WHERE geom IS NOT NULL {month_filter}
        ),
        binned AS (
            SELECT
                ROUND(lon::numeric / 0.05) * 0.05 AS lon_bin,
                ROUND(lat::numeric / 0.05) * 0.05 AS lat_bin,
                COUNT(*) AS cnt
            FROM pts
            WHERE lon BETWEEN 29 AND 42 AND lat BETWEEN -7 AND 5
            GROUP BY lon_bin, lat_bin
        )
        SELECT lon_bin AS lon, lat_bin AS lat, cnt
        FROM binned
        ORDER BY cnt DESC
        LIMIT 8000
    """
    try:
        cur.execute(query)
        rows = [dict(r) for r in cur.fetchall()]
        conn.close()

        if not rows:
            return {"data": [], "breaks": []}

        counts = [int(r["cnt"]) for r in rows]
        nonzero = sorted(c for c in counts if c > 0)
        n = len(nonzero)

        if n >= 3:
            p33 = nonzero[int(n * 0.33)]
            p66 = nonzero[int(n * 0.66)]
        else:
            p33 = nonzero[0] if nonzero else 1
            p66 = nonzero[-1] if nonzero else 2

        breaks = [0, p33, p66]

        result = []
        for r in rows:
            c = int(r["cnt"])
            if c == 0:
                cls = 1
            elif c <= p33:
                cls = 2
            elif c <= p66:
                cls = 3
            else:
                cls = 4
            result.append({"lat": float(r["lat"]), "lon": float(r["lon"]), "count": c, "class": cls})

        return {"data": result, "breaks": breaks, "season": season}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        try:
            conn.close()
        except:
            pass


# ─── SOC LAYER ─────────────────────────────────────────────────────────────────
@app.get("/api/soc/sample")
def soc_sample():
    conn = get_conn()
    cur = dict_cursor(conn)
    query = """
        SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS soc_value
        FROM (
            SELECT (ST_PixelAsCentroids(rast)).*
            FROM o_4_soc_dataset
        ) AS pixels
        WHERE val > 0 AND val < 500
        LIMIT 4000
    """
    try:
        cur.execute(query)
        rows = cur.fetchall()
        return {"data": [dict(r) for r in rows]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ─── RAINFALL LAYER ────────────────────────────────────────────────────────────
@app.get("/api/rainfall/sample")
def rainfall_sample(season: Optional[str] = "janmar"):
    if season not in SEASONS:
        raise HTTPException(status_code=400, detail="Invalid season")
    table = f"{SEASONS[season]['table_suffix']}Rainfall"
    conn = get_conn()
    cur = dict_cursor(conn)
    query = f"""
        SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS rainfall_mm
        FROM (
            SELECT (ST_PixelAsCentroids(rast)).*
            FROM "{table}"
        ) AS pixels
        WHERE val > 0
        LIMIT 3000
    """
    try:
        cur.execute(query)
        return {"data": [dict(r) for r in cur.fetchall()]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ─── NDVI LAYER (season-specific) ──────────────────────────────────────────────
@app.get("/api/ndvi/sample")
def ndvi_sample(season: Optional[str] = "janmar"):
    # Try seasonal table, fall back to janmar
    table = NDVI_TABLES.get(season, "ndvi_janmar_2019")
    conn = get_conn()
    cur = dict_cursor(conn)
    for tbl in [table, "ndvi_janmar_2019"]:
        try:
            query = f"""
                SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS ndvi
                FROM (
                    SELECT (ST_PixelAsCentroids(rast)).*
                    FROM {tbl}
                ) AS pixels
                WHERE val BETWEEN -1 AND 1
                LIMIT 3000
            """
            cur.execute(query)
            rows = cur.fetchall()
            conn.close()
            return {"data": [dict(r) for r in rows], "table": tbl}
        except:
            continue
    conn.close()
    return {"data": [], "table": None}


# ─── DEM LAYER ─────────────────────────────────────────────────────────────────
@app.get("/api/dem/sample")
def dem_sample():
    conn = get_conn()
    cur = dict_cursor(conn)
    query = """
        SELECT ST_X(geom) AS lon, ST_Y(geom) AS lat, val AS elevation
        FROM (
            SELECT (ST_PixelAsCentroids(rast)).*
            FROM (SELECT rast FROM "30mdem" LIMIT 6) tiles
        ) AS pixels
        WHERE val > 0
        LIMIT 2000
    """
    try:
        cur.execute(query)
        return {"data": [dict(r) for r in cur.fetchall()]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()


# ─── SEASONAL ANALYSIS (dissertation Fig 6, 13, 14, 15) ───────────────────────
@app.get("/api/analysis/seasonal")
def seasonal_analysis():
    """
    For each season: GPS count, avg SOC at migration points, avg NDVI at grazed cells,
    rainfall from dissertation Table 2, and lagged NDVI (T+1 season NDVI at T grazed locations).
    """
    conn = get_conn()
    cur = dict_cursor(conn)
    results = []

    for season_key, season_info in SEASONS.items():
        months = season_info["months"]
        mf = f"EXTRACT(MONTH FROM TO_DATE(date, 'DD/MM/YYYY')) = ANY(ARRAY{months}::int[])"

        # GPS count
        try:
            cur.execute(f"SELECT COUNT(*) AS n FROM migration_data WHERE geom IS NOT NULL AND {mf}")
            n = int(cur.fetchone()["n"])
        except:
            n = 0

        avg_soc = avg_ndvi = avg_elev = avg_ndvi_lag = None

        try:
            # Sample up to 300 migration points for this season
            cur.execute(f"""
                WITH sample AS (
                    SELECT {MIGRATION_TO_WGS84} AS gw
                    FROM migration_data
                    WHERE geom IS NOT NULL AND {mf}
                    LIMIT 300
                )
                SELECT
                    AVG((SELECT ST_Value(s.rast, p.gw) FROM o_4_soc_dataset s
                         WHERE ST_Intersects(s.rast, p.gw) LIMIT 1)) AS avg_soc,
                    AVG((SELECT ST_Value(d.rast, p.gw) FROM "30mdem" d
                         WHERE ST_Intersects(d.rast, p.gw) LIMIT 1)) AS avg_elev
                FROM sample p
            """)
            row = cur.fetchone()
            if row:
                avg_soc = float(row["avg_soc"]) if row["avg_soc"] else None
                avg_elev = float(row["avg_elev"]) if row["avg_elev"] else None
        except Exception as e:
            print(f"SOC/elev error ({season_key}): {e}")

        # NDVI at grazed locations for this season (current NDVI)
        ndvi_tbl = NDVI_TABLES.get(season_key, "ndvi_janmar_2019")
        try:
            cur.execute(f"""
                WITH sample AS (
                    SELECT {MIGRATION_TO_WGS84} AS gw
                    FROM migration_data
                    WHERE geom IS NOT NULL AND {mf}
                    LIMIT 200
                )
                SELECT AVG((SELECT ST_Value(n.rast, p.gw) FROM {ndvi_tbl} n
                     WHERE ST_Intersects(n.rast, p.gw) LIMIT 1)) AS avg_ndvi
                FROM sample p
            """)
            r2 = cur.fetchone()
            avg_ndvi = float(r2["avg_ndvi"]) if r2 and r2["avg_ndvi"] else None
        except Exception as e:
            print(f"NDVI error ({season_key}): {e}")

        # Lagged NDVI: T+1 season's NDVI at T's grazed locations
        next_s = NEXT_SEASON[season_key]
        ndvi_tbl_lag = NDVI_TABLES.get(next_s, "ndvi_janmar_2019")
        try:
            cur.execute(f"""
                WITH sample AS (
                    SELECT {MIGRATION_TO_WGS84} AS gw
                    FROM migration_data
                    WHERE geom IS NOT NULL AND {mf}
                    LIMIT 200
                )
                SELECT AVG((SELECT ST_Value(n.rast, p.gw) FROM {ndvi_tbl_lag} n
                     WHERE ST_Intersects(n.rast, p.gw) LIMIT 1)) AS avg_ndvi_lag
                FROM sample p
            """)
            r3 = cur.fetchone()
            avg_ndvi_lag = float(r3["avg_ndvi_lag"]) if r3 and r3["avg_ndvi_lag"] else None
        except Exception as e:
            print(f"Lagged NDVI error ({season_key}): {e}")

        results.append({
            "season": season_key,
            "label": season_info["label"],
            "n_points": n,
            "avg_soc": round(avg_soc, 2) if avg_soc else None,
            "avg_rainfall": SEASON_RAINFALL[season_key],
            "avg_elevation": round(avg_elev, 1) if avg_elev else None,
            "avg_ndvi": round(avg_ndvi, 4) if avg_ndvi else None,
            "avg_ndvi_lag": round(avg_ndvi_lag, 4) if avg_ndvi_lag else None,
            "next_season_label": SEASONS[next_s]["label"],
        })

    conn.close()
    return {"data": results}


# ─── GRAZED vs UNGRAZED SOC (dissertation Figs 19, 20, 22, 23, Table 5, 6) ───
@app.get("/api/analysis/grazed-vs-ungrazed")
def grazed_vs_ungrazed():
    """
    Core dissertation analysis:
    1. Classify migration data cells into GRAZED (class 2-4) or UNGRAZED (class 1)
       using full multi-year dataset (2013-2018, all months)
    2. For each class: avg SOC, avg elevation
    3. Compute residual SOC = SOC - (slope*elev + intercept) after OLS regression
    4. Returns per-point data for scatter + summary stats for Table 6
    """
    conn = get_conn()
    cur = dict_cursor(conn)

    # Step 1: Get all binned cells with grazing count and raster values
    query = f"""
        WITH all_pts AS (
            SELECT
                ROUND(ST_X({MIGRATION_TO_WGS84})::numeric / 0.1) * 0.1 AS lon_bin,
                ROUND(ST_Y({MIGRATION_TO_WGS84})::numeric / 0.1) * 0.1 AS lat_bin,
                COUNT(*) AS cnt
            FROM migration_data
            WHERE geom IS NOT NULL
            GROUP BY lon_bin, lat_bin
        ),
        grid AS (
            SELECT
                lon_bin, lat_bin, cnt,
                ST_SetSRID(ST_MakePoint(lon_bin, lat_bin), 4326) AS pt
            FROM all_pts
            WHERE lon_bin BETWEEN 29 AND 42 AND lat_bin BETWEEN -7 AND 5
        )
        SELECT
            g.lon_bin, g.lat_bin, g.cnt,
            (SELECT ST_Value(s.rast, g.pt) FROM o_4_soc_dataset s
             WHERE ST_Intersects(s.rast, g.pt) LIMIT 1) AS soc,
            (SELECT ST_Value(d.rast, g.pt) FROM "30mdem" d
             WHERE ST_Intersects(d.rast, g.pt) LIMIT 1) AS elevation
        FROM grid g
        ORDER BY g.cnt DESC
        LIMIT 2000
    """
    try:
        cur.execute(query)
        rows = [dict(r) for r in cur.fetchall()]
        conn.close()
    except Exception as e:
        conn.close()
        raise HTTPException(status_code=500, detail=str(e))

    if not rows:
        return {"data": [], "summary": {}}

    # Step 2: Jenks-like classification using percentiles on non-zero counts
    counts = [int(r["cnt"]) for r in rows]
    nonzero = sorted(c for c in counts if c > 0)
    n = len(nonzero)
    p33 = nonzero[int(n*0.33)] if n >= 3 else 1
    p66 = nonzero[int(n*0.66)] if n >= 3 else 2

    valid = []
    for r in rows:
        c = int(r["cnt"])
        soc = float(r["soc"]) if r["soc"] and 0 < float(r["soc"]) < 500 else None
        elev = float(r["elevation"]) if r["elevation"] else None
        if soc is None or elev is None:
            continue
        cls = 1 if c == 0 else (2 if c <= p33 else (3 if c <= p66 else 4))
        grazed = cls > 1
        valid.append({
            "lon": float(r["lon_bin"]),
            "lat": float(r["lat_bin"]),
            "count": c,
            "class": cls,
            "grazed": grazed,
            "soc": soc,
            "elevation": elev,
        })

    if not valid:
        return {"data": [], "summary": {}}

    # Step 3: OLS regression SOC ~ Elevation to get residuals
    xs = [v["elevation"] for v in valid]
    ys = [v["soc"] for v in valid]
    slope, intercept, r2 = linear_regression(xs, ys)

    for v in valid:
        if slope is not None:
            predicted = slope * v["elevation"] + intercept
            v["residual_soc"] = round(v["soc"] - predicted, 4)
        else:
            v["residual_soc"] = 0.0

    # Step 4: Summary stats per group
    def group_stats(pts):
        if not pts:
            return {}
        socs = [p["soc"] for p in pts]
        resids = [p["residual_soc"] for p in pts]
        elevs = [p["elevation"] for p in pts]
        n = len(pts)
        return {
            "n": n,
            "mean_soc": round(sum(socs)/n, 2),
            "mean_residual_soc": round(sum(resids)/n, 2),
            "mean_elevation": round(sum(elevs)/n, 1),
            "sd_soc": round(math.sqrt(sum((x-sum(socs)/n)**2 for x in socs)/n), 2) if n > 1 else 0,
            "sd_residual_soc": round(math.sqrt(sum((x-sum(resids)/n)**2 for x in resids)/n), 2) if n > 1 else 0,
        }

    grazed_pts = [v for v in valid if v["grazed"]]
    ungrazed_pts = [v for v in valid if not v["grazed"]]

    # By class
    class_stats = {}
    for cls in [1, 2, 3, 4]:
        label_map = {1: "Ungrazed", 2: "Mild", 3: "Moderate", 4: "Heavy"}
        cls_pts = [v for v in valid if v["class"] == cls]
        class_stats[label_map[cls]] = group_stats(cls_pts)

    # Pearson correlation: residual SOC ~ grazing class
    rsc = pearson([v["count"] for v in valid], [v["residual_soc"] for v in valid])
    rraw = pearson([v["count"] for v in valid], [v["soc"] for v in valid])

    return {
        "data": valid,
        "summary": {
            "grazed": group_stats(grazed_pts),
            "ungrazed": group_stats(ungrazed_pts),
            "by_class": class_stats,
            "ols_slope": slope,
            "ols_intercept": intercept,
            "ols_r2": r2,
            "corr_count_residual_soc": rsc,
            "corr_count_raw_soc": rraw,
            "breaks": [0, p33, p66],
        }
    }


# ─── NDVI GRAZED vs UNGRAZED (dissertation Figs 14, 15, 17) ──────────────────
@app.get("/api/analysis/ndvi-grazing")
def ndvi_grazing():
    """
    For each season: mean NDVI at grazed vs ungrazed cells.
    Also returns the lagged comparison (grazing in T → NDVI in T+1).
    Matches dissertation Figs 14, 15, 17.
    """
    conn = get_conn()
    cur = dict_cursor(conn)

    # Get long-term grazed/ungrazed classification from full dataset
    try:
        cur.execute(f"""
            WITH all_pts AS (
                SELECT
                    ROUND(ST_X({MIGRATION_TO_WGS84})::numeric / 0.1) * 0.1 AS lon_bin,
                    ROUND(ST_Y({MIGRATION_TO_WGS84})::numeric / 0.1) * 0.1 AS lat_bin,
                    COUNT(*) AS cnt
                FROM migration_data WHERE geom IS NOT NULL
                GROUP BY lon_bin, lat_bin
            )
            SELECT lon_bin, lat_bin, cnt,
                ST_SetSRID(ST_MakePoint(lon_bin, lat_bin), 4326) AS pt
            FROM all_pts
            WHERE lon_bin BETWEEN 29 AND 42 AND lat_bin BETWEEN -7 AND 5
        """)
        cells = [dict(r) for r in cur.fetchall()]
    except Exception as e:
        conn.close()
        raise HTTPException(status_code=500, detail=str(e))

    counts = [int(c["cnt"]) for c in cells]
    nonzero = sorted(x for x in counts if x > 0)
    n = len(nonzero)
    p33 = nonzero[int(n*0.33)] if n >= 3 else 1

    # Classify grazed (cnt > p33) vs ungrazed (cnt <= p33 or 0)
    grazed_pts = [c for c in cells if int(c["cnt"]) > p33]
    ungrazed_pts = [c for c in cells if int(c["cnt"]) <= p33]

    results = []
    for season_key in SEASONS:
        ndvi_tbl = NDVI_TABLES.get(season_key, "ndvi_janmar_2019")
        lag_tbl = NDVI_TABLES.get(NEXT_SEASON[season_key], "ndvi_janmar_2019")

        def get_ndvi_for_pts(pts, tbl, limit=150):
            if not pts:
                return None
            # Sample a subset for performance
            sample = pts[:limit]
            pt_list = ",".join(f"ST_SetSRID(ST_MakePoint({p['lon_bin']},{p['lat_bin']}),4326)" for p in sample)
            try:
                cur.execute(f"""
                    WITH pts AS (SELECT unnest(ARRAY[{pt_list}]) AS gw)
                    SELECT AVG((SELECT ST_Value(n.rast, p.gw) FROM {tbl} n
                                WHERE ST_Intersects(n.rast, p.gw) LIMIT 1)) AS avg_ndvi
                    FROM pts p
                """)
                row = cur.fetchone()
                return float(row["avg_ndvi"]) if row and row["avg_ndvi"] else None
            except:
                return None

        grazed_ndvi = get_ndvi_for_pts(grazed_pts, ndvi_tbl)
        ungrazed_ndvi = get_ndvi_for_pts(ungrazed_pts, ndvi_tbl)
        grazed_lag = get_ndvi_for_pts(grazed_pts, lag_tbl)
        ungrazed_lag = get_ndvi_for_pts(ungrazed_pts, lag_tbl)

        results.append({
            "season": season_key,
            "label": SEASONS[season_key]["label"],
            "next_label": SEASONS[NEXT_SEASON[season_key]]["label"],
            "rainfall": SEASON_RAINFALL[season_key],
            "grazed_ndvi": round(grazed_ndvi, 4) if grazed_ndvi else None,
            "ungrazed_ndvi": round(ungrazed_ndvi, 4) if ungrazed_ndvi else None,
            "grazed_ndvi_lag": round(grazed_lag, 4) if grazed_lag else None,
            "ungrazed_ndvi_lag": round(ungrazed_lag, 4) if ungrazed_lag else None,
        })

    conn.close()
    return {"data": results}


# ─── CORRELATION DATA (raw scatter points) ────────────────────────────────────
@app.get("/api/analysis/correlation")
def correlation_data():
    conn = get_conn()
    cur = dict_cursor(conn)
    query = f"""
        WITH grid AS (
            SELECT
                ROUND(ST_X({MIGRATION_TO_WGS84})::numeric, 1) AS lon_bin,
                ROUND(ST_Y({MIGRATION_TO_WGS84})::numeric, 1) AS lat_bin,
                COUNT(*) AS grazing_count,
                ST_SetSRID(ST_MakePoint(
                    AVG(ST_X({MIGRATION_TO_WGS84})),
                    AVG(ST_Y({MIGRATION_TO_WGS84}))
                ), 4326) AS center_pt
            FROM migration_data WHERE geom IS NOT NULL
            GROUP BY lon_bin, lat_bin
        )
        SELECT
            lon_bin, lat_bin, grazing_count,
            (SELECT ST_Value(s.rast, g.center_pt) FROM o_4_soc_dataset s
             WHERE ST_Intersects(s.rast, g.center_pt) LIMIT 1) AS soc,
            (SELECT ST_Value(r.rast, g.center_pt) FROM janmar2019Rainfall r
             WHERE ST_Intersects(r.rast, g.center_pt) LIMIT 1) AS rainfall,
            (SELECT ST_Value(d.rast, g.center_pt) FROM "30mdem" d
             WHERE ST_Intersects(d.rast, g.center_pt) LIMIT 1) AS elevation,
            (SELECT ST_Value(n.rast, g.center_pt) FROM ndvi_janmar_2019 n
             WHERE ST_Intersects(n.rast, g.center_pt) LIMIT 1) AS ndvi
        FROM grid g
        WHERE lon_bin BETWEEN 29 AND 42 AND lat_bin BETWEEN -7 AND 5
        ORDER BY grazing_count DESC
        LIMIT 800
    """
    try:
        cur.execute(query)
        return {"data": [dict(r) for r in cur.fetchall()]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        conn.close()
