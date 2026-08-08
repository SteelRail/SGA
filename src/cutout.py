"""Brick coadd on disk -> SCI/IVAR/MASK sample with WCS and PSF metadata.

A sample is one multi-extension FITS, identical for galaxy and background
frames so the consumer never learns which driver produced a file:

    hdul['SCI']   (3, H, W) float32, grz, nanomaggies, WCS in the header
    hdul['IVAR']  (3, H, W) float32, 1/nanomaggy²
    hdul['MASK']  (H, W) int16, DR9 MASKBITS

The header carries PSF_G/R/Z (FWHM in arcsec read from the psfsize planes
at the sample's own position, so the consumer's forward model takes its
kernel width from the data rather than a config constant), PIXSCALE,
BRICK, and provenance identifiers supplied by the driver (SGA_ID for
galaxies; the parent SGA_ID, offset and position angle for backgrounds).

No flux rescaling, no stretch, no clipping: values stay in nanomaggies
and negative sky pixels survive. Stamps that run past the brick edge are
zero-filled with IVAR 0 and the NPRIMARY bit set, so the invalid region
is marked rather than the frame rejected; how much of a frame must be
valid is the drivers' single structural criterion.
"""

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

PIXSCALE = 0.262  # arcsec per pixel, DR9 brick grid
BANDS = ("g", "r", "z")


def brick_dir(root, brickname):
    """Mirror directory of one brick's files under the data root."""
    return Path(root) / "bricks" / brickname[:3] / brickname


# Fill values for the off-brick region of edge-clipped stamps.
FILL_SCI = 0.0
FILL_IVAR = 0.0
FILL_MASK = 1  # MASKBITS bit 0, NPRIMARY


class Coadd:
    """One brick's coadd planes in memory, ready to be cut many times."""

    def __init__(self, directory, brickname):
        self.brickname = brickname
        self.directory = Path(directory)

        def plane(kind, band=None):
            suffix = f"{kind}-{band}" if band else kind
            path = self.directory / f"legacysurvey-{brickname}-{suffix}.fits.fz"
            with fits.open(path) as hdul:
                hdu = hdul[1] if hdul[0].data is None else hdul[0]
                return hdu.data, hdu.header

        images, headers = zip(*(plane("image", band) for band in BANDS), strict=True)
        self.image = np.stack(images).astype(np.float32)
        self.invvar = np.stack(
            [plane("invvar", band)[0] for band in BANDS]
        ).astype(np.float32)
        self.psfsize = np.stack(
            [plane("psfsize", band)[0] for band in BANDS]
        ).astype(np.float32)
        maskbits, _ = plane("maskbits")
        self.maskbits = maskbits.astype(np.int16)
        self.wcs = WCS(headers[0])
        self.shape = self.image.shape[1:]

    def contains(self, ra_deg, dec_deg, half_px=0):
        """True when the stamp box around the position touches the pixel grid at all."""
        x, y = self.wcs.world_to_pixel_values(ra_deg, dec_deg)
        ny, nx = self.shape
        return (-half_px <= x < nx + half_px) and (-half_px <= y < ny + half_px)

    def psf_fwhm(self, ra_deg, dec_deg, size_px):
        """Per-band PSF FWHM (arcsec) at a position, from the psfsize planes.

        The value at the centre pixel; where that carries no coverage
        (psfsize 0), the median over the covered part of the stamp box.
        """
        x, y = self.wcs.world_to_pixel_values(ra_deg, dec_deg)
        col = int(np.clip(round(float(x)), 0, self.shape[1] - 1))
        row = int(np.clip(round(float(y)), 0, self.shape[0] - 1))
        half = size_px // 2
        rows = slice(max(row - half, 0), min(row + half + 1, self.shape[0]))
        cols = slice(max(col - half, 0), min(col + half + 1, self.shape[1]))
        fwhm = []
        for plane in self.psfsize:
            value = float(plane[row, col])
            if not value > 0:
                box = plane[rows, cols]
                covered = box[box > 0]
                value = float(np.median(covered)) if covered.size else 0.0
            fwhm.append(value)
        return fwhm


def cut(brick, ra_deg, dec_deg, size_px):
    """Cut one (SCI, IVAR, MASK) stamp of `size_px` pixels centred on a position.

    Returns a dict with the three planes, the stamp WCS (carried by
    Cutout2D so a sky position round-trips to the same pixel), the
    per-band PSF FWHM, and the fraction of pixels with coverage in all
    three bands (`valid_frac`, the structural quality measure).
    """
    position = brick.wcs.world_to_pixel_values(ra_deg, dec_deg)
    position = (float(position[0]), float(position[1]))
    size = (int(size_px), int(size_px))

    def plane2d(data, fill):
        stamp = Cutout2D(data, position, size, wcs=brick.wcs,
                         mode="partial", fill_value=fill)
        return stamp

    sci = [plane2d(brick.image[i], FILL_SCI) for i in range(len(BANDS))]
    ivar = [plane2d(brick.invvar[i], FILL_IVAR).data for i in range(len(BANDS))]
    mask = plane2d(brick.maskbits, FILL_MASK)

    sci_stack = np.stack([s.data for s in sci]).astype(np.float32)
    ivar_stack = np.stack(ivar).astype(np.float32)
    mask_plane = np.asarray(mask.data, dtype=np.int16)

    return {
        "sci": sci_stack,
        "ivar": ivar_stack,
        "mask": mask_plane,
        "wcs": sci[0].wcs,
        "psf_fwhm": brick.psf_fwhm(ra_deg, dec_deg, size_px),
        "valid_frac": float(np.mean(np.all(ivar_stack > 0, axis=0))),
        "ra": float(ra_deg),
        "dec": float(dec_deg),
    }


def write_sample(path, sample, brickname, provenance, layers=None):
    """Write one sample to `path` in the fixed contract.

    `provenance` is a dict of extra header cards — (value, comment)
    tuples — identifying the sample's origin. `layers`, when given, is an
    (H, W) uint8 bit-packed plane appended as a fourth extension LAYERS
    (bit meanings in its header); the first three extensions are the
    contract and do not depend on it.
    """
    header = fits.Header()
    header.update(sample["wcs"].to_header())
    header["EXTNAME"] = "SCI"
    header["BUNIT"] = ("nanomaggy", "AB zeropoint 22.5")
    header["PIXSCALE"] = (PIXSCALE, "arcsec per pixel")
    header["BRICK"] = (brickname, "DR9 brick")
    header["BANDS"] = (",".join(BANDS), "planes along axis 3")
    header["RA0"] = (sample["ra"], "requested centre RA (deg)")
    header["DEC0"] = (sample["dec"], "requested centre Dec (deg)")
    for band, fwhm in zip(BANDS, sample["psf_fwhm"], strict=True):
        header[f"PSF_{band.upper()}"] = (fwhm, "PSF FWHM at this position (arcsec)")
    header["VALIDFRC"] = (sample["valid_frac"], "fraction covered in all bands")
    for key, value in provenance.items():
        header[key] = value

    ivar_header = fits.Header()
    ivar_header["EXTNAME"] = "IVAR"
    ivar_header["BUNIT"] = "1/nanomaggy^2"
    mask_header = fits.Header()
    mask_header["EXTNAME"] = "MASK"
    mask_header["COMMENT"] = "DR9 MASKBITS; legacysurvey.org/dr9/bitmasks"

    hdus = [
        fits.PrimaryHDU(data=sample["sci"], header=header),
        fits.ImageHDU(data=sample["ivar"], header=ivar_header),
        fits.ImageHDU(data=sample["mask"], header=mask_header),
    ]
    if layers is not None:
        layer_header = fits.Header()
        layer_header["EXTNAME"] = "LAYERS"
        layer_header["BIT0"] = ("INVALID", "no coverage in some band")
        layer_header["BIT1"] = ("BRIGHT", "bright object per MASKBITS")
        layer_header["BIT2"] = ("SOURCE", "tractor detection footprint")
        layer_header["BIT3"] = ("GALAXY", "SGA ellipse at D26")
        hdus.append(fits.ImageHDU(data=layers, header=layer_header))

    path = Path(path)
    path.parent.mkdir(exist_ok=True, parents=True)
    fits.HDUList(hdus).writeto(path, overwrite=True)
