# Gaia Sidecar Generator

Create human-readable CSV sidecar files with Gaia bright-star information for your galaxy cutout images.

## Overview

`make_gaia_sidecars.py` reads a folder of JPEG cutouts and generates an associated CSV file for each image containing:

- **Source IDs** from Gaia DR3
- **Coordinates**: RA, Dec (degrees)
- **Magnitude**: G-band (mag)
- **Separation**: Angular distance and pixel distance from cutout center
- **Pixel position**: x, y coordinates within the cutout

## Usage

### Basic Example

```bash
python make_gaia_sidecars.py \
  --cutout_folder ./galaxy_cutouts \
  --gaia_catalog ./bright_stars_gaia_dr3_b5.fits
```

This creates, for each cutout like:
```
sga_1234_RA150.1234_Dec2.3456_size90.0.jpeg
```

A corresponding sidecar:
```
sga_1234_RA150.1234_Dec2.3456_size90.0_gaia.csv
```

### With Custom Settings

```bash
python make_gaia_sidecars.py \
  --cutout_folder ./data/SGA \
  --gaia_catalog ./bright_stars_gaia_dr3_b5.fits \
  --pixscale 0.262 \
  --gaia_maglim 16.5 \
  --overwrite
```

## Command-Line Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--cutout_folder` | Yes | - | Folder containing JPEG cutouts |
| `--gaia_catalog` | Yes | - | Path to Gaia bright-star FITS catalog |
| `--pixscale` | No | 0.262 | Pixel scale (arcsec/pixel); must match your download settings |
| `--gaia_maglim` | No | 17.0 | G-band magnitude limit for stars in sidecars |
| `--overwrite` | No | False | Overwrite existing sidecar files |

## Output Format

### Sidecar CSV Columns

```
source_id,ra_deg,dec_deg,Gmag,sep_arcsec,sep_pix,x_pix,y_pix
4358647842309261312,150.1245,2.3401,14.23,12.4,47.4,224.3,151.2
4358647842309261440,150.1289,2.3567,15.67,8.1,31.0,186.5,203.8
4358647842309261568,150.0998,2.3302,16.12,14.7,56.1,342.1,180.4
```

### Column Definitions

- **source_id**: Gaia DR3 unique identifier (if present in input catalog)
- **ra_deg**: Right ascension in degrees
- **dec_deg**: Declination in degrees
- **Gmag**: Gaia G-band magnitude
- **sep_arcsec**: Angular separation from cutout center (arcsec)
- **sep_pix**: Separation in pixels
- **x_pix**: X pixel coordinate (0 = left edge, increases to the right/west)
- **y_pix**: Y pixel coordinate (0 = bottom edge, increases upward/north)

## Supported Filename Formats

The script automatically parses both:

**Galaxy mode:**
```
sga_2_RA228.3771_Dec5.4232_size44.5.jpeg
```

**Fields mode:**
```
sga_2_field_0.0_RA228.4599_Dec5.4232_size44.5.jpeg
```

The parser extracts RA, Dec, and size_arcsec from the filename.

## Pixel Coordinate System

- **Origin**: (0, 0) at the lower-left corner of the cutout
- **X axis**: Increases to the right (pointing west on the sky)
- **Y axis**: Increases upward (pointing north on the sky)
- **Valid range**: 0 ≤ x_pix < size_px, 0 ≤ y_pix < size_px

Pixel coordinates are computed using small-angle approximation (tangent plane projection) centered on the cutout center coordinates.

## Workflow Example

### Step 1: Download galaxy cutouts

```bash
python load_sga.py \
  --mode galaxy \
  --n_limit 1000 \
  --download_folder ./galaxies
```

### Step 2: Create Gaia bright-star catalog

```bash
python create_bright_stars_table.py \
  --gal-lat-min 5.0 \
  --output-fits bright_stars.fits
```

### Step 3: Generate sidecars

```bash
python make_gaia_sidecars.py \
  --cutout_folder ./galaxies \
  --gaia_catalog bright_stars.fits \
  --gaia_maglim 16.0
```

Each galaxy cutout now has a CSV sidecar listing nearby Gaia stars!

## Why Sidecars?

**Advantages of sidecar files:**

1. **Decoupled**: One file per cutout, easy to process independently
2. **Human-readable**: Plain-text CSV, can be inspected in any text editor or spreadsheet
3. **Lightweight**: ~1 KB per file vs. re-querying the database
4. **Reproducible**: Exact same information for long-term archival
5. **Scalable**: Parallel processing of sidecars without database load

## Performance

- **Speed**: ~100–500 cutouts per second (depends on star density)
- **I/O**: Reads FITS catalog once, writes one CSV per cutout with stars
- **Memory**: Full Gaia catalog loaded into RAM (~5–10 GB for bright stars)

For 100,000 cutouts with typical star density, expect ~3–10 minutes total runtime.

## Troubleshooting

### "Could not parse RA/Dec/size from filename"

Your JPEG filename doesn't match the expected format. Ensure it contains:
```
_RA<number>_Dec<number>_size<number>.jpeg
```

### "No Gaia stars found in any cutout"

Check:
1. Is `--gaia_maglim` too bright? (e.g., 17.0 is typical)
2. Does your Gaia catalog actually contain stars? Load with:
   ```python
   from astropy.table import Table
   gaia = Table.read('bright_stars.fits')
   print(len(gaia), gaia.colnames)
   ```

### Sidecar files are empty

Stars may be outside the cutout boundaries due to pixscale mismatch. Verify:
```bash
python load_sga.py --help | grep pixscale
python make_gaia_sidecars.py --help | grep pixscale
```
Must match!

## Integration with Your Pipeline

The sidecar CSVs are designed to work seamlessly with:
- **ML/CV pipelines**: Load CSVs alongside JPEGs for star detection training
- **Morphological analysis**: Mask out bright stars before galaxy analysis
- **Astrometric validation**: Cross-check cutout coordinates against Gaia
- **Quality control**: Identify crowded fields before detailed analysis

## References

- Gaia DR3: [ESA Gaia Release](https://www.cosmos.esa.int/web/gaia/dr3)
- Legacy Survey: [DR9 Images](https://www.legacysurvey.org/dr9/)
