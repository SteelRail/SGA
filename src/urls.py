"""Remote data layout: every NERSC static-server file the pipeline fetches."""

DR9_BASE = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9"
SGA_BASE = "https://portal.nersc.gov/project/cosmo/data/sga/2020"

CATALOGUES = {
    "SGA-2020.fits": f"{SGA_BASE}/SGA-2020.fits",
    "survey-bricks-dr9-north.fits.gz":
        f"{DR9_BASE}/north/survey-bricks-dr9-north.fits.gz",
    "survey-bricks-dr9-south.fits.gz":
        f"{DR9_BASE}/south/survey-bricks-dr9-south.fits.gz",
}


def brick_files(brickname, hemisphere):
    """URLs of one brick's coadd files and tractor catalogue, keyed by filename."""
    coadd = f"{DR9_BASE}/{hemisphere}/coadd/{brickname[:3]}/{brickname}"
    urls = {}
    for band in ("g", "r", "z"):
        for kind in ("image", "invvar", "psfsize"):
            name = f"legacysurvey-{brickname}-{kind}-{band}.fits.fz"
            urls[name] = f"{coadd}/{name}"
    name = f"legacysurvey-{brickname}-maskbits.fits.fz"
    urls[name] = f"{coadd}/{name}"
    name = f"tractor-{brickname}.fits"
    urls[name] = f"{DR9_BASE}/{hemisphere}/tractor/{brickname[:3]}/{name}"
    return urls


def brick_checksums(brickname, hemisphere):
    """URL of the sha256 manifest published next to one brick's coadds."""
    coadd = f"{DR9_BASE}/{hemisphere}/coadd/{brickname[:3]}/{brickname}"
    return (f"{coadd}/legacysurvey_dr9_{hemisphere}_coadd_"
            f"{brickname[:3]}_{brickname}.sha256sum")


def psf_stamps(group_ra, group_name):
    """URLs of one SGA-2020 group's pixelised PSF stamps, keyed by filename."""
    base = f"{SGA_BASE}/data/{int(group_ra):03d}/{group_name}"
    return {
        f"{group_name}-largegalaxy-psf-{band}.fits.fz":
            f"{base}/{group_name}-largegalaxy-psf-{band}.fits.fz"
        for band in ("g", "r", "z")
    }
