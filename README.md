# SGA training-data builder

Research tools for downloading DESI Legacy Survey DR9 imaging of SGA-2020
galaxies — calibrated grz cutouts with inverse variance, masks and PSF
metadata — as training data for GRAF.

Both classes of sample — galaxies and the backgrounds they are trained
against — are cut from the same L1 brick coadds, fetched from the static
NERSC file servers (never the viewer's cutout renderer, never the SGA
mosaics). Frames are not rejected for containing bright objects; four
independent mask layers mark pixels instead, so no criterion separates
the classes on anything but the galaxy. Values stay in nanomaggies with
negative sky intact, and each frame's header carries the PSF FWHM read
from the survey's `psfsize` maps at that position.

## Setup

```bash
uv sync                                    # .venv with pinned dependencies
uv run python scripts/fetch_catalogues.py  # reference catalogues (~0.8 GB, one-off)
```

Everything large — brick mirror, samples, PSF stamps — lands under
`data/sga` (every script takes `--root` to point elsewhere); nothing
under `data/` is tracked.

## Usage

```bash
# 1. Galaxy samples: select -> fetch bricks -> cut (each stage resumable;
#    start with a spread subset, drop --sample for the full selection)
uv run python scripts/fetch_galaxies.py --sample 200

# 2. Background samples for the same selection, with survival report
uv run python scripts/fetch_backgrounds.py --sample 200

# 3. One-off PSF calibration set (SGA pixelised PSFs -> Moffat profile)
uv run python scripts/fetch_psf_library.py --n-groups 300
```

Each sample is one FITS file, identical for both drivers:

```python
hdul = fits.open('.../sga_1179597.fits')
hdul['SCI']     # (3, H, W) float32, grz, nanomaggies, WCS in header
hdul['IVAR']    # (3, H, W) float32, 1/nanomaggy^2
hdul['MASK']    # (H, W) int16, DR9 MASKBITS
hdul['LAYERS']  # (H, W) uint8: invalid / bright / source / galaxy bits
# header: PSF_G/R/Z (arcsec), PIXSCALE, BRICK, SGA_ID (or PARENT+OFFSET)
```

Measured on a 200-galaxy spread: ~4 MB per sample (~0.4 TB for the full
z < 0.05 selection), ~90 MB per brick (~6 TB to mirror all 69k bricks;
bricks are discardable once samples are cut). Background survival is
~85%, with ~86% of accepted backgrounds sharing the parent's brick.

## Layout

| path | purpose |
|---|---|
| `src/urls.py` | remote data layout: NERSC file URLs |
| `src/geometry.py` | sky geometry: ellipse radii, tangent-plane offsets |
| `src/bricks.py` | survey-bricks tables: footprint test, position -> brick |
| `src/catalog.py` | SGA-2020: selection, ellipses, full-atlas clearance |
| `src/select.py` | galaxy target selection composing catalogue and bricks |
| `src/fetch.py` | static-HTTP mirroring: session reuse, resume, checksums |
| `src/cutout.py` | brick coadd -> SCI/IVAR/MASK sample with WCS and PSF |
| `src/mask.py` | the four mask layers, per-sample verdict record |
| `scripts/fetch_catalogues.py` | one-off reference catalogue download |
| `scripts/fetch_galaxies.py` | SGA positions -> bricks -> samples |
| `scripts/fetch_backgrounds.py` | stepped positions -> bricks -> samples |
| `scripts/fetch_psf_library.py` | SGA PSF stamps -> parametric profile fit |
| `data/catalogues/` | SGA-2020 and survey-bricks tables (untracked) |
| `notebooks/01_inspect_samples.ipynb` | tour of the sample contract |
| `notebooks/02_background_purity.ipynb` | accepted/rejected stamps by eye, survival |
| `legacy/` | superseded viewer-cutout scripts, kept for reference |

Selection criteria, separations and mask definitions are documented in
the module and script docstrings, next to the parameters that set them.
Linting is manual: `uv run ruff check src scripts`.
