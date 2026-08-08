# Download Pipeline Manager (`run_download_pipeline.py`)

## Overview

`run_download_pipeline.py` is a comprehensive manager script that orchestrates the complete workflow for downloading and processing galaxy and field cutout images from DECALS, with automatic integrity checking and Gaia bright-star sidecar generation.

## Pipeline Stages

The pipeline executes **5 sequential stages** with automatic error handling:

### Stage 1: Download Galaxy Central Cutouts
- Downloads SGA galaxy images centered on each target
- Uses DECALS DR9 footprint filtering to skip out-of-survey galaxies
- Default: adaptive image sizing based on galaxy semi-major axis
- Output: `data/cutouts/galaxies/*.jpeg`

### Stage 2: Galaxy Integrity Check & Retry
- Checks each downloaded galaxy cutout for validity
- Identifies broken files (< 500 bytes, typically HTML errors from failed API calls)
- Re-downloads all broken and missing files
- Continues to **Stage 3** only after all galaxy files are verified

### Stage 3: Download Background Field Cutouts
- For each galaxy, downloads N background field positions at large radius (~20× D26)
- Used as control/background samples
- Offset calculation prevents overlaps with galaxy positions
- Output: `data/cutouts/fields/*.jpeg`

### Stage 4: Field Integrity Check & Retry
- Same integrity checking as Stage 2 but for field cutouts
- Catches network failures, API rate-limit issues, etc.
- Re-downloads broken/missing field files

### Stage 5: Generate Gaia Bright-Star Sidecars
- Runs separately on `galaxies/` and `fields/` folders
- For each JPEG cutout, finds all Gaia bright stars (g < 17.0 mag) within the frame
- Generates `*_gaia.csv` sidecar with:
  - Star positions (RA, Dec, pixel coordinates)
  - Gaia G-band magnitudes
  - Angular separation from cutout center
- Output: `data/cutouts/galaxies/*_gaia.csv` and `data/cutouts/fields/*_gaia.csv`

## Usage

### Basic Run (100 galaxies, 7 fields each)
```bash
python run_download_pipeline.py --n_galaxies 100
```

### Custom Configuration
```bash
python run_download_pipeline.py \
  --n_galaxies 500 \
  --n_fields 5 \
  --cutouts_folder ./my_cutouts \
  --num_workers 8 \
  --pixscale 0.262 \
  --size -1
```

### Selective Re-runs
Skip early stages if data already downloaded:
```bash
# Skip downloads, only run retry + sidecars
python run_download_pipeline.py --n_galaxies 100 --skip_galaxies --skip_fields

# Only redo galaxy retry and sidecars
python run_download_pipeline.py --n_galaxies 100 --skip_fields

# Only generate sidecars (if cutouts already validated)
python run_download_pipeline.py --n_galaxies 100 --skip_galaxies --skip_fields
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--n_galaxies` | 100 | Number of SGA galaxies to download |
| `--n_fields` | 7 | Background field positions per galaxy |
| `--cutouts_folder` | `./data/cutouts` | Parent folder for galaxy/ and fields/ subfolders |
| `--catalog` | `./data/catalogues/SGA-2020.fits` | SGA galaxy catalog path |
| `--gaia_catalog` | `./data/catalogues/gaia_bright_stars_g15_b5.fits` | Gaia bright-star catalog (g<15, \|b\|≥5°) |
| `--pixscale` | 0.262 | Pixel scale (arcsec/pixel) - must match DECALS |
| `--size` | -1 | Image size: -1=adaptive, or fixed pixels |
| `--num_workers` | 4 | Parallel download workers (3-4 recommended) |
| `--bricks_file` | `./data/catalogues/survey-bricks-dr9-south.fits` | DECALS bricks (southern hemisphere) |
| `--bricks_file_north` | `./data/catalogues/survey-bricks-dr9-north.fits` | DECALS bricks (northern hemisphere) |
| `--broken_size_threshold` | 500 | Min file size (bytes) to consider valid |
| `--skip_galaxies` | False | Skip galaxy download and retry stages |
| `--skip_fields` | False | Skip field download and retry stages |
| `--skip_sidecars` | False | Skip Gaia sidecar generation |

## Output Structure

```
data/cutouts/
├── galaxies/
│   ├── sga_1_RA123.45_Dec45.67_size120.5.jpeg
│   ├── sga_1_RA123.45_Dec45.67_size120.5_gaia.csv
│   ├── sga_2_RA...jpeg
│   └── ...
└── fields/
    ├── sga_1_field_1_RA180.12_Dec10.34_size100.0.jpeg
    ├── sga_1_field_1_RA180.12_Dec10.34_size100.0_gaia.csv
    ├── sga_1_field_2_RA...jpeg
    └── ...
```

### Gaia Sidecar Format (CSV)
```
source_id,ra_deg,dec_deg,Gmag,sep_arcsec,sep_pix,x_pix,y_pix
5897382957289374,123.456,45.678,12.3,45.2,172.6,512.1,256.8
5897382957289375,123.467,45.689,13.1,48.5,185.3,498.7,243.2
```

Columns:
- `source_id` - Gaia DR3 unique identifier
- `ra_deg`, `dec_deg` - Star's celestial coordinates
- `Gmag` - Gaia G-band magnitude
- `sep_arcsec` - Angular separation from cutout center
- `sep_pix` - Separation in pixels
- `x_pix`, `y_pix` - Star position within cutout frame (0,0 = bottom-left)

## Test Results

Pipeline tested successfully on **100 SGA galaxies**:

```
Configuration:
  N galaxies:        100
  N fields/galaxy:   7
  Pixel scale:       0.262 arcsec/pixel
  Workers:           4

Results:
  Galaxy cutouts:    100 JPEG files, 6 sidecars
  Field cutouts:     719 JPEG files, 67 sidecars
  Total:             819 JPEG files, 73 sidecars
  Execution time:    ~10 minutes
```

**Why only 6 galaxy sidecars?**
- Gaia bright stars (g<15, |b|≥5) are rare in most DECALS fields
- Most of the 100 galaxy cutouts have NO stars falling within the frame boundaries
- 67 field sidecars generated because field positions were randomly placed across sky

## Catalog Files

### Gaia Bright Stars Catalog
- **File:** `gaia_bright_stars_g15_b5.fits`
- **Location:** `~/data/catalogues/` (shared with collaborators)
- **Symlink:** `./data/catalogues/gaia_bright_stars_g15_b5.fits` → central copy
- **Size:** 916 MB
- **Filters applied:**
  - G magnitude < 15.0
  - Galactic latitude |b| ≥ 5.0°
- **Generation:** Run `python create_bright_stars_table.py --output-fits ./data/catalogues/gaia_bright_stars_g15_b5.fits`
- **Rows:** 26,682,475 bright stars

### DECALS Footprint Files
- **South:** `survey-bricks-dr9-south.fits` (73 MB)
- **North:** `survey-bricks-dr9-north.fits` (27 MB)
- **Location:** `./data/catalogues/`
- **Purpose:** Define valid regions for cutout downloads (prevents wasted API calls to out-of-survey coordinates)

### SGA Catalog
- **File:** `SGA-2020.fits` (symlink to external location)
- **Size:** 683 MB
- **Rows:** 104,906 galaxies
- **Source:** Legacy Survey galaxy catalog (Simons Observatory)

## Key Features

1. **Automatic Failure Recovery**
   - Identifies broken downloads (< 500 byte threshold = HTML error)
   - Re-downloads with exponential backoff
   - Won't proceed past any stage until all files validated

2. **Parallel Downloads**
   - Configurable worker pool (default 4, max 8)
   - Respects API rate limits
   - Adaptive speed based on network response

3. **DECALS Footprint Filtering**
   - Pre-filters galaxies/fields against survey bricks
   - Skips ~6% of out-of-survey targets automatically
   - Reduces failed API calls and improves efficiency

4. **Progress Reporting**
   - Real-time progress bars for each stage
   - Step-by-step output with command echoing
   - Summary statistics at completion

5. **Flexible Configuration**
   - Separate galaxy/field folders for easy organization
   - Skip individual stages for iterative workflows
   - Custom image sizes, pixel scales, worker counts

6. **Reproducible Outputs**
   - Consistent naming convention (RA, Dec, size in filename)
   - Gaia sidecars with precise pixel coordinates
   - Full pipeline summary in terminal output

## Typical Workflow

```bash
# 1. Download and validate 500 galaxies + 7 fields each
python run_download_pipeline.py --n_galaxies 500 --num_workers 4

# 2. Later: redo only sidecars if Gaia catalog updated
python run_download_pipeline.py --n_galaxies 500 --skip_galaxies --skip_fields

# 3. Download different field set (10 fields instead of 7)
python run_download_pipeline.py --n_galaxies 500 --n_fields 10 --skip_galaxies
```

## Troubleshooting

**Issue:** "Gaia catalog not found"
- **Solution:** Run `python create_bright_stars_table.py` to generate it

**Issue:** "DECALS bricks file not found"
- **Solution:** Check paths in `./data/catalogues/` or update via `--bricks_file` argument

**Issue:** Most cutouts have 0 sidecars
- **Solution:** Normal if few bright stars in those sky regions. Check Gaia catalog filters (g<15, |b|≥5°)

**Issue:** High failure rate during downloads
- **Solution:** Reduce `--num_workers` (try 2-3 instead of 4), increase retry threshold with `--broken_size_threshold`

**Issue:** Out of memory during field downloads
- **Solution:** Reduce `--n_galaxies` or `--n_fields` to lower total files

## See Also

- `load_sga.py` - Underlying download script (supports galaxy/fields/retry modes)
- `make_gaia_sidecars.py` - Gaia sidecar generation (called internally by pipeline)
- `create_bright_stars_table.py` - Gaia catalog query tool
- README.md - Project overview and setup instructions
