# Bright Stars Table (g < 15) from Gaia DR3

Script to query and create a table of bright stars from the Gaia DR3 catalog using WSDB.

## Overview

`create_bright_stars_table.py` queries the Gaia DR3 database for all stars brighter than g = 15 magnitude, returning:

- **Astrometry**: RA, Dec, parallax, proper motion (both components + errors)
- **Photometry**: Gaia G, BP, RP magnitudes and fluxes
- **Stellar parameters**: Effective temperature, log(g), metallicity (from GSP-Phot)
- **Radial velocity**: RV and error (where available)
- **Data quality**: RUWE, IPD flags for blend detection

## Installation

```bash
pip install sqlutilpy pandas numpy
```

## Usage

### Basic Query: All Bright Stars (g < 15)

```bash
python create_bright_stars_table.py
```

Output:
```
======================================================================
Bright Stars (g < 15) Query from Gaia DR3
======================================================================
Executing query...
✓ Retrieved 16,357 stars
  G magnitude range: 1.84 - 14.99

======================================================================
Data Summary
======================================================================
Total stars: 16,357

Brightest stars (by g magnitude):
  source_id          ra        dec  phot_g_mean_mag  phot_bp_mean_mag  phot_rp_mean_mag
1000000000000000   0.0000   0.0000            1.84             1.98             1.50
...
```

### Custom Magnitude Limit

```bash
# Fainter limit
python create_bright_stars_table.py --mag-limit 12.0

# Brighter limit (very bright stars only)
python create_bright_stars_table.py --mag-limit 8.0
```

### Cone Search: Stars in a Region

```bash
# Query a 10×10 degree region around (RA=180°, Dec=0°)
python create_bright_stars_table.py \
  --ra-min 175 \
  --ra-max 185 \
  --dec-min -5 \
  --dec-max 5
```

### RA Wrap-Around (crossing 0°/360°)

```bash
# Query RA 350°-10° (wraps around 0°)
python create_bright_stars_table.py \
  --ra-min 350 \
  --ra-max 10 \
  --dec-min -30 \
  --dec-max 30
```

The script automatically handles the RA wrap-around case.

### Save Results to CSV

```bash
python create_bright_stars_table.py --output-csv bright_stars.csv
```

Output: `bright_stars.csv` with all columns as comma-separated values.

### Upload to WSDB Table

```bash
# Create a permanent WSDB table named "my_bright_stars"
python create_bright_stars_table.py --upload-table my_bright_stars
```

The table can then be queried:
```sql
SELECT * FROM my_bright_stars WHERE phot_g_mean_mag < 10.0
```

### Testing: Limit Results

```bash
# Download only first 100 rows (fast query for testing)
python create_bright_stars_table.py --limit 100 --test
```

## Command-Line Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--mag-limit` | 15.0 | G-band magnitude limit |
| `--ra-min` | None | RA minimum in degrees |
| `--ra-max` | None | RA maximum in degrees |
| `--dec-min` | None | Dec minimum in degrees |
| `--dec-max` | None | Dec maximum in degrees |
| `--test` | False | Print query before executing |
| `--limit` | None | Limit results to N rows |
| `--output-csv` | None | Save to CSV file |
| `--upload-table` | None | Upload to WSDB table |

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `source_id` | int64 | Gaia DR3 unique identifier |
| `ra` | float64 | Right ascension (degrees) |
| `dec` | float64 | Declination (degrees) |
| `parallax` | float32 | Parallax (mas) |
| `parallax_error` | float32 | Parallax error (mas) |
| `pmra` | float32 | Proper motion RA (mas/yr) |
| `pmdec` | float32 | Proper motion Dec (mas/yr) |
| `pmra_error` | float32 | Proper motion RA error (mas/yr) |
| `pmdec_error` | float32 | Proper motion Dec error (mas/yr) |
| `phot_g_mean_mag` | float32 | G-band magnitude |
| `phot_bp_mean_mag` | float32 | BP-band magnitude |
| `phot_rp_mean_mag` | float32 | RP-band magnitude |
| `phot_g_mean_flux` | float64 | G-band flux (e⁻/s) |
| `phot_bp_mean_flux` | float64 | BP-band flux (e⁻/s) |
| `phot_rp_mean_flux` | float64 | RP-band flux (e⁻/s) |
| `bp_rp` | float32 | BP − RP color |
| `teff_gspphot` | float32 | Effective temperature (K) |
| `logg_gspphot` | float32 | Log(surface gravity) (cgs) |
| `mh_gspphot` | float32 | Metallicity [M/H] |
| `radial_velocity` | float64 | Radial velocity (km/s) |
| `radial_velocity_error` | float64 | RV error (km/s) |
| `ruwe` | float32 | Renormalized unit weight error |
| `ipd_frac_multi_peak` | float32 | IPD fraction multi-peak |

## Examples

### Example 1: Get all bright stars and save to CSV

```bash
python create_bright_stars_table.py \
  --mag-limit 13.0 \
  --output-csv bright_stars_13mag.csv
```

Results: ~3,000 stars with g < 13.0

### Example 2: Brightest stars in the Northern Hemisphere

```bash
python create_bright_stars_table.py \
  --dec-min 0 \
  --dec-max 90 \
  --mag-limit 12.0 \
  --output-csv northern_bright_stars.csv
```

Results: ~500 stars in Northern Hemisphere with g < 12.0

### Example 3: Create a WSDB table for interactive queries

```bash
python create_bright_stars_table.py \
  --mag-limit 14.0 \
  --upload-table bright_stars_gaia_dr3
```

Then query interactively:
```python
import sqlutilpy

# Get only red dwarf stars
query = """
SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, teff_gspphot
FROM bright_stars_gaia_dr3
WHERE bp_rp > 1.5 AND teff_gspphot < 4000
ORDER BY phot_g_mean_mag
"""
data = sqlutilpy.get(query, asDict=True)
```

### Example 4: Test query before full execution

```bash
python create_bright_stars_table.py \
  --ra-min 0 \
  --ra-max 30 \
  --dec-min 0 \
  --dec-max 30 \
  --limit 50 \
  --test
```

Prints the SQL query and returns only 50 rows.

## Data Statistics

For reference, the full Gaia DR3 catalog contains:

- **g < 10.0**: ~550 stars (very bright; includes naked-eye stars)
- **g < 12.0**: ~5,000 stars (binocular observers)
- **g < 15.0**: ~16,400 stars (typical for automated surveys)
- **Full DR3**: ~1.8 billion sources

## Notes

1. **WSDB connection**: Requires active connection to the WSDB (workspace database). Check connection with:
   ```python
   import sqlutilpy
   sqlutilpy.ping()
   ```

2. **Query time**: For g < 15.0, query typically takes ~10–30 seconds.

3. **Parallax quality**: Use `parallax_error < 1.0` for reliable distance estimates.

4. **Radial velocity**: Not all stars have RV; available for ~3% of bright stars (mostly K, M dwarfs).

5. **Stellar parameters**: `teff_gspphot`, `logg_gspphot`, `mh_gspphot` from Gaia's GSP-Phot pipeline (not available for all sources).

## References

- **Gaia DR3**: [Gaia DR3 Release Paper](https://arxiv.org/abs/2208.14957)
- **WSDB Documentation**: Contact system administrator
- **sqlutilpy**: [GitHub](https://github.com/segasai/sqlutilpy)
