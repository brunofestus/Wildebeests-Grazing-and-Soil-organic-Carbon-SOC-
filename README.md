# Wildebeest Migration & SOC Analysis — Full-Stack App

## Project Structure
```
wildebeest-app/
├── backend/
│   ├── main.py           # FastAPI backend
│   ├── requirements.txt  # Python deps
│   └── .env.example      # DB config template
└── frontend/
    └── index.html        # Single-file frontend (open in browser)
```

---

## Setup & Run

### 1. Configure Database

```bash
cd backend
cp .env.example .env
```

Edit `.env` with your PostgreSQL credentials:
```
DB_HOST=localhost
DB_PORT=5432
DB_NAME=your_database_name
DB_USER=postgres
DB_PASSWORD=your_password
```

### 2. Install Python Dependencies

```bash
cd backend
pip install -r requirements.txt
```

Or with a virtual environment (recommended):
```bash
cd backend
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

### 3. Start the Backend

```bash
cd backend
uvicorn main:app --reload --port 8000
```

You should see:
```
INFO: Uvicorn running on http://127.0.0.1:8000
```

Test it: open http://localhost:8000/api/health in your browser — should return `{"status":"ok"}`

### 4. Open the Frontend

Simply open `frontend/index.html` in your browser:
- Double-click the file, OR
- Serve it: `python -m http.server 3000` from the `frontend/` folder, then go to http://localhost:3000

---

## API Endpoints

| Endpoint | Description |
|---|---|
| `GET /api/health` | DB connection check |
| `GET /api/migration/density?season=janmar` | Migration density grid |
| `GET /api/soc/sample` | SOC raster sampled to points |
| `GET /api/dem/sample` | DEM raster (resampled) |
| `GET /api/rainfall/sample?season=janmar` | Seasonal rainfall |
| `GET /api/ndvi/sample` | NDVI (Jan–Mar 2019) |
| `GET /api/analysis/seasonal` | Season-by-season aggregations |
| `GET /api/analysis/correlation` | Grid-cell level data for scatter plots |
| `GET /api/analysis/regression` | Pearson correlations & summary stats |
| `GET /api/study-area` | Study area boundary GeoJSON |

### Season keys
`janmar` · `aprmay` · `junjuly` · `augoct` · `novdec`

---

## Troubleshooting

**"Backend offline" shown in the app**
→ Make sure uvicorn is running on port 8000

**Empty map layers**
→ Check the browser console (F12) for errors; the API endpoint `/api/migration/density` should return data

**Raster queries are slow**
→ This is normal for large PostGIS rasters; consider creating raster overviews:
```sql
SELECT AddRasterConstraints('public', 'o_4_soc_dataset', 'rast');
SELECT ST_CreateOverview('o_4_soc_dataset'::regclass, 'rast', 4);
```

**"relation does not exist" errors**
→ Make sure table names match exactly (case-sensitive). Update the `SEASONS` dict in `main.py` if needed.
