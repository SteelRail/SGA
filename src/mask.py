"""The four mask layers and the per-sample verdict record.

Masks mark pixels; they never reject frames. A bright star in a frame is
masked and kept, because a rejection criterion applied to one class and
not the other separates the classes on something other than the galaxy.
The only structural criterion is coverage (`valid_frac` against
DEFAULT_COVERAGE_MIN), and both drivers must apply the same threshold —
which is why it lives here and not in either driver.

Each layer is computed independently and stored as one bit of the LAYERS
plane, so any criterion can be revised without recomputing the others:

    bit 0  INVALID  no coverage in at least one band (invvar == 0; also
                    covers chip gaps and the off-brick fill of
                    edge-clipped stamps)
    bit 1  BRIGHT   damaged pixels per MASKBITS: SATUR_G/R/Z and
                    ALLMASK_G/R/Z. The DR9 BRIGHT/MEDIUM bits are
                    deliberately excluded — they are magnitude-radius
                    proximity circles around catalogued stars, not
                    pixel damage, and swallow whole frames near bright
                    stars; a consumer that wants them reads the raw
                    MASK plane
    bit 2  SOURCE   tractor detection footprints from TYPE and
                    SHAPE_R/E1/E2 (the frame's own SGA galaxies excluded
                    — ref_cat L3 — so the galaxy layer alone owns them)
    bit 3  GALAXY   SGA-2020 ellipses at D26, every catalogued galaxy
                    overlapping the frame

Source footprints: extended types (REX/EXP/DEV/SER) are masked to
`radius_factor` times the half-light radius SHAPE_R along the fitted
ellipse; point types (PSF/DUP) to `psf_factor` times the frame's r-band
PSF FWHM. Tractor ellipticity components are e1 = e cos(2b),
e2 = e sin(2b) with b the position angle North towards East, verified
against the SGA catalogue PAs of the tractor L3 rows.

MASKBITS bit numbers follow legacysurvey.org/dr9/bitmasks.
"""

import numpy as np

from .cutout import PIXSCALE
from .geometry import elliptical_radius

# legacysurvey.org/dr9/bitmasks
MASKBITS = {
    "NPRIMARY": 0,
    "BRIGHT": 1,
    "SATUR_G": 2,
    "SATUR_R": 3,
    "SATUR_Z": 4,
    "ALLMASK_G": 5,
    "ALLMASK_R": 6,
    "ALLMASK_Z": 7,
    "WISEM1": 8,
    "WISEM2": 9,
    "BAILOUT": 10,
    "MEDIUM": 11,
    "GALAXY": 12,
    "CLUSTER": 13,
}

BRIGHT_BITS = sum(
    1 << MASKBITS[name]
    for name in ("SATUR_G", "SATUR_R", "SATUR_Z",
                 "ALLMASK_G", "ALLMASK_R", "ALLMASK_Z")
)
GALAXY_BIT = 1 << MASKBITS["GALAXY"]

LAYER_INVALID = 1
LAYER_BRIGHT = 2
LAYER_SOURCE = 4
LAYER_GALAXY = 8

# The single structural threshold, shared by both drivers: a frame whose
# all-band coverage falls below this is unusable and is not written.
DEFAULT_COVERAGE_MIN = 0.8


def invalid_layer(ivar):
    """Pixels without trustworthy data: no coverage in at least one band."""
    return ~np.all(ivar > 0, axis=0)


def bright_layer(mask_plane):
    """Pixels with actual damage per MASKBITS: saturation and ALLMASK."""
    return (mask_plane & BRIGHT_BITS) != 0


def _sky_jacobian(wcs, x0, y0):
    """2x2 matrix mapping pixel steps to (east, north) offsets in arcsec."""
    ra = np.array([x0, x0 + 1.0, x0])
    dec = np.array([y0, y0, y0 + 1.0])
    world_ra, world_dec = wcs.pixel_to_world_values(ra, dec)
    cos_dec = np.cos(np.radians(world_dec[0]))
    d_east = ((world_ra[1:] - world_ra[0] + 180.0) % 360.0 - 180.0) * cos_dec * 3600.0
    d_north = (world_dec[1:] - world_dec[0]) * 3600.0
    return np.array([[d_east[0], d_east[1]], [d_north[0], d_north[1]]])


def _paint_ellipse(layer, x0, y0, a_arcsec, ba, pa_deg, jacobian):
    """Set True inside the ellipse (semi-major `a_arcsec`) centred at pixel (x0, y0)."""
    height, width = layer.shape
    reach = int(np.ceil(a_arcsec / PIXSCALE)) + 1
    x_lo, x_hi = int(np.floor(x0)) - reach, int(np.ceil(x0)) + reach + 1
    y_lo, y_hi = int(np.floor(y0)) - reach, int(np.ceil(y0)) + reach + 1
    x_lo, x_hi = max(x_lo, 0), min(x_hi, width)
    y_lo, y_hi = max(y_lo, 0), min(y_hi, height)
    if x_lo >= x_hi or y_lo >= y_hi:
        return
    ys, xs = np.mgrid[y_lo:y_hi, x_lo:x_hi]
    d_east = jacobian[0, 0] * (xs - x0) + jacobian[0, 1] * (ys - y0)
    d_north = jacobian[1, 0] * (xs - x0) + jacobian[1, 1] * (ys - y0)
    inside = elliptical_radius(d_east, d_north, pa_deg, ba) <= a_arcsec
    layer[y_lo:y_hi, x_lo:x_hi] |= inside


def source_layer(shape, wcs, tractor, psf_fwhm_r,
                 radius_factor=2.0, psf_factor=1.25, min_radius=1.0):
    """Per-detection footprints of every tractor source except SGA galaxies."""
    layer = np.zeros(shape, dtype=bool)
    # plain-array copies: fancy-indexing the FITS record array per sample
    # costs ~30 ms; these copies cost well under 1 ms
    keep = np.asarray(tractor["ref_cat"]) != "L3"
    if not keep.any():
        return layer
    rows = {column: np.asarray(tractor[column])[keep] for column in
            ("ra", "dec", "type", "shape_r", "shape_e1", "shape_e2")}
    x, y = wcs.world_to_pixel_values(rows["ra"], rows["dec"])

    types = np.char.strip(rows["type"].astype(str))
    point = np.isin(types, ("PSF", "DUP"))
    shape_r = np.array(rows["shape_r"], dtype=float)
    radius = np.where(point,
                      np.maximum(psf_factor * psf_fwhm_r, min_radius),
                      np.maximum(radius_factor * shape_r, min_radius))

    e1 = np.array(rows["shape_e1"], dtype=float)
    e2 = np.array(rows["shape_e2"], dtype=float)
    e = np.hypot(e1, e2)
    ba = np.where(point | (e <= 0), 1.0, (1.0 - e) / (1.0 + e))
    pa = 0.5 * np.degrees(np.arctan2(e2, e1))

    height, width = shape
    margin = radius / PIXSCALE
    near = (x > -margin) & (x < width + margin) & (y > -margin) & (y < height + margin)

    jacobian = _sky_jacobian(wcs, (width - 1) / 2.0, (height - 1) / 2.0)
    for i in np.flatnonzero(near):
        _paint_ellipse(layer, float(x[i]), float(y[i]), float(radius[i]),
                       float(ba[i]), float(pa[i]), jacobian)
    return layer


def galaxy_layer(shape, wcs, catalog, ra0, dec0, margin=1.0):
    """SGA-2020 ellipses at `margin` * r26, for every galaxy reaching the frame."""
    layer = np.zeros(shape, dtype=bool)
    height, width = shape
    half_diag = 0.5 * np.hypot(height, width) * PIXSCALE
    indices = catalog.overlapping(ra0, dec0, half_diag, margin=margin)
    if not len(indices):
        return layer
    jacobian = _sky_jacobian(wcs, (width - 1) / 2.0, (height - 1) / 2.0)
    for index in indices:
        r26, ba, pa = catalog.ellipse(int(index))
        x, y = wcs.world_to_pixel_values(catalog._ra[index], catalog._dec[index])
        _paint_ellipse(layer, float(x), float(y), margin * r26, ba, pa, jacobian)
    return layer


def compute_layers(sample, catalog, tractor,
                   radius_factor=2.0, psf_factor=1.25, galaxy_margin=1.0):
    """All four layers for one cut sample; returns (layers_uint8, verdict).

    The verdict records each layer's pixel fraction separately, plus the
    fraction of DR9's own GALAXY maskbit, so selection criteria can be
    revisited from the manifest without touching pixels.
    """
    shape = sample["mask"].shape
    wcs = sample["wcs"]
    layers = np.zeros(shape, dtype=np.uint8)

    invalid = invalid_layer(sample["ivar"])
    bright = bright_layer(sample["mask"])
    source = source_layer(shape, wcs, tractor, sample["psf_fwhm"][1],
                          radius_factor=radius_factor, psf_factor=psf_factor)
    galaxy = galaxy_layer(shape, wcs, catalog, sample["ra"], sample["dec"],
                          margin=galaxy_margin)

    layers |= np.uint8(LAYER_INVALID) * invalid
    layers |= np.uint8(LAYER_BRIGHT) * bright
    layers |= np.uint8(LAYER_SOURCE) * source
    layers |= np.uint8(LAYER_GALAXY) * galaxy

    verdict = {
        "valid_frac": sample["valid_frac"],
        "bright_frac": float(bright.mean()),
        "source_frac": float(source.mean()),
        "galaxy_frac": float(galaxy.mean()),
        "galaxy_bit_frac": float(((sample["mask"] & GALAXY_BIT) != 0).mean()),
    }
    return layers, verdict
