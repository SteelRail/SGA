# load_sga.py

Download JPEG cutouts from DESI Legacy Survey DR9 for SGA-2020 galaxies.

## Overview

`load_sga.py` downloads optical image cutouts from the Legacy Survey DR9 for nearby galaxies in the SGA-2020 catalog. Supports two download modes:

1. **Galaxy mode**: Central galaxy cutouts
2. **Fields mode**: Background field cutouts around each galaxy

## Data Source

**SGA-2020 (Siena Galaxy Atlas)**
- 383,620 nearby galaxies (z < 0.05 = 104,906 galaxies)
- Reference: Moustakas et al. 2023 ([arXiv:2307.04888](https://arxiv.org/abs/2307.04888))
- Catalog: https://www.legacysurvey.org/sga/sga2020/

**Imaging**: DESI Legacy Survey DR9 optical (grz bands)

## Installation

```bash
source ~/Work/venvs/.venv/bin/activate
pip install astropy requests numpy tqdm
```

## Usage

### Galaxy Mode: Download Central Galaxy Cutouts

Download JPEG images of galaxy centers with adaptive sizing based on galaxy diameter (D26):

```bash
python load_sga.py \
  --catalog /Users/vasilybelokurov/data/catalogues/SGA-2020.fits \
  --mode galaxy \
  --download_folder ./galaxy_cutouts \
  --num_workers 8 \
  --n_limit 1000
```

Output filenames:
```
sga_2_RA228.3771_Dec5.4232_size44.5.jpeg
sga_24_RA247.4245_Dec40.2482_size48.9.jpeg
...
```

### Fields Mode: Download Background Field Cutouts

Download background field cutouts around each galaxy. Field positions are uniformly distributed on a circle centered on the galaxy:

```bash
python load_sga.py \
  --catalog /Users/vasilybelokurov/data/catalogues/SGA-2020.fits \
  --mode fields \
  --n_away 7 \
  --na_away 20 \
  --download_folder ./field_cutouts \
  --num_workers 8 \
  --n_limit 100
```

Output filenames (one per field position, per galaxy):
```
sga_2_field_0.0_RA228.4599_Dec5.4232_size44.5.jpeg      (angle 0°)
sga_2_field_51.4_RA228.4206_Dec5.5120_size44.5.jpeg     (angle 51.4°)
sga_2_field_102.9_RA228.3162_Dec5.6765_size44.5.jpeg    (angle 102.9°)
...
sga_24_field_0.0_RA247.5431_Dec40.2482_size48.9.jpeg
...
```

## Command-Line Arguments

### Common Arguments

- `--catalog PATH` (default: `./data/catalogues/SGA-2020.fits`)
  - Path to SGA-2020 FITS catalog

- `--download_folder PATH` (default: `./data/SGA`)
  - Output directory for downloaded images

- `--zcut FLOAT` (default: `0.05`)
  - Redshift cutoff for galaxy selection

- `--num_workers INT` (default: `16`)
  - Number of parallel download workers

- `--chunksize INT` (default: `32`)
  - Chunk size for parallel processing

- `--n_limit INT` (default: `None`)
  - Limit to N galaxies (for testing)

- `--layer STR` (default: `ls-dr9`)
  - Legacy Survey layer for cutouts

- `--pixscale FLOAT` (default: `0.262`)
  - Pixel scale in arcseconds/pixel

### Sizing Arguments

- `--size INT` (default: `-1`)
  - Cutout size in pixels:
    - `-1`: Adaptive sizing (recommended)
    - `N > 0`: Fixed size (all cutouts N×N pixels)
  - Adaptive formula: `size = na × (D26/2 in arcsec) / pixscale`
  - `na=3.0` by default (3× the galaxy semi-major axis)

- `--na FLOAT` (default: `3.0`)
  - Multiplier for adaptive sizing (galaxy mode)
  - Larger values = larger cutouts

### Mode-Specific Arguments

#### Galaxy Mode (default)

```bash
python load_sga.py --mode galaxy [common args] [sizing args]
```

No mode-specific arguments; uses adaptive sizing by default.

#### Fields Mode

```bash
python load_sga.py --mode fields [common args] [sizing args] --n_away N --na_away NA
```

- `--n_away INT` (default: `7`)
  - Number of field positions per galaxy
  - Positions are uniformly distributed on a circle

- `--na_away FLOAT` (default: `20.0`)
  - Field offset as multiplier of galaxy radius
  - Field radius = `na_away × D26/2` (in arcsec)
  - Field cutout size = galaxy cutout size (adaptive)

## Mathematical Details

### Galaxy Mode Sizing

For each galaxy with diameter D26 (in arcminutes):

```
size_arcsec = na × (D26 / 2) × 60  [convert arcmin → arcsec, diameter → radius]
size_pixels = size_arcsec / pixscale
```

Default: `na = 3.0`, `pixscale = 0.262` arcsec/pixel

Example (D26 = 4 arcmin):
```
size_arcsec = 3.0 × 2 × 60 = 360 arcsec
size_pixels = 360 / 0.262 ≈ 1374 pixels
```

### Fields Mode Positioning

For each galaxy at (RA_gal, Dec_gal) with diameter D26:

1. **Field circle radius** (in degrees):
   ```
   radius_arcsec = na_away × (D26 / 2) × 60
   radius_degrees = radius_arcsec / 3600
   ```

2. **Field positions** (uniformly distributed):
   - For i = 0, 1, ..., n_away-1:
   ```
   angle_i = 360 × i / n_away  [degrees]
   RA_field = RA_gal + (radius_deg / cos(Dec_gal)) × cos(angle_rad)
   Dec_field = Dec_gal + radius_deg × sin(angle_rad)
   ```
   - Declination is accounted for when computing RA offset
   - Declination clamped to [-90°, +90°]

3. **Field cutout size**: Same as galaxy cutout size (adaptive with `na = 3.0`)

Example (D26 = 4 arcmin, na_away = 20, n_away = 7):
```
radius_arcsec = 20 × 2 × 60 = 2400 arcsec ≈ 0.67°
Field positions at angles: 0°, 51.4°, 102.9°, 154.3°, 205.7°, 257.1°, 308.6°
```

## Output Filename Format

### Galaxy Mode

```
sga_<SGA_ID>_RA<RA>_Dec<DEC>_size<ARCSEC>.jpeg
```

- `SGA_ID`: Galaxy ID from SGA-2020 catalog (integer)
- `RA`: Right ascension in degrees (4 decimal places)
- `DEC`: Declination in degrees (4 decimal places)
- `ARCSEC`: Cutout size in arcseconds (1 decimal place)

Example: `sga_2_RA228.3771_Dec5.4232_size44.5.jpeg`

### Fields Mode

```
sga_<SGA_ID>_field_<ANGLE>_RA<RA>_Dec<DEC>_size<ARCSEC>.jpeg
```

- `SGA_ID`: Parent galaxy ID
- `ANGLE`: Position angle on field circle in degrees (1 decimal place)
- `RA`, `DEC`, `ARCSEC`: As above

Example: `sga_2_field_0.0_RA228.4599_Dec5.4232_size44.5.jpeg`

## Examples

### Test: Download 10 galaxies

```bash
python load_sga.py --n_limit 10 --mode galaxy --num_workers 4
```

### Production: Download all galaxies with 8 workers

```bash
python load_sga.py --catalog /path/to/SGA-2020.fits \
  --mode galaxy \
  --download_folder ./galaxy_data \
  --num_workers 8 \
  --zcut 0.05
```

Expected: ~104,906 galaxies, ~0.5 TB data

### Background fields: 5 positions per galaxy for 500 galaxies

```bash
python load_sga.py --n_limit 500 \
  --mode fields \
  --n_away 5 \
  --na_away 15 \
  --download_folder ./field_data \
  --num_workers 8
```

Expected: 500 × 5 = 2,500 field cutouts
