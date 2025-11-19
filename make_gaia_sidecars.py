#!/usr/bin/env python3
"""
Create Gaia bright-star sidecar CSV files for galaxy cutouts.

For each JPEG cutout in a folder, this script:
  1. Parses RA, Dec, and size_arcsec from the filename
  2. Queries a Gaia bright-star catalog for stars inside the cutout
  3. Computes pixel coordinates and angular separations
  4. Writes a human-readable CSV sidecar with star information

Sidecar files are saved as <image_basename>_gaia.csv in the same folder.
"""

import os
import re
import argparse

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create Gaia bright-star sidecar CSV files for galaxy cutouts."
    )
    parser.add_argument(
        "--cutout_folder",
        type=str,
        required=True,
        help="Folder containing JPEG cutouts.",
    )
    parser.add_argument(
        "--gaia_catalog",
        type=str,
        required=True,
        help="Path to Gaia bright-star catalog (FITS with ra, dec, phot_g_mean_mag).",
    )
    parser.add_argument(
        "--pixscale",
        type=float,
        default=0.262,
        help="Pixel scale in arcsec per pixel (default: 0.262).",
    )
    parser.add_argument(
        "--gaia_maglim",
        type=float,
        default=17.0,
        help="G-band magnitude limit for Gaia stars in sidecars (default: 17.0).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing _gaia.csv sidecars.",
    )
    return parser.parse_args()


def load_gaia_catalog(path):
    """Load Gaia bright-star catalog and construct SkyCoord array.

    Expected columns: ra, dec, phot_g_mean_mag.
    Optionally: source_id, gal_l, gal_b.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Gaia catalog file not found: {path}")

    gaia = Table.read(path)
    required_cols = ["ra", "dec", "phot_g_mean_mag"]
    if not all(col in gaia.colnames for col in required_cols):
        raise ValueError(
            "Gaia catalog must contain columns: ra, dec, phot_g_mean_mag"
        )

    gaia_coord = SkyCoord(gaia["ra"] * u.deg, gaia["dec"] * u.deg)
    print(f"Loaded Gaia catalog with {len(gaia)} entries from {path}")
    return gaia, gaia_coord


def parse_cutout_filename(filename):
    """Parse RA, Dec, and size_arcsec from a cutout filename.

    Works for both:
      sga_{id}_RA{ra:.4f}_Dec{dec:.4f}_size{size_arcsec:.1f}.jpeg
      sga_{id}_field_{angle:.1f}_RA{ra:.4f}_Dec{dec:.4f}_size{size_arcsec:.1f}.jpeg

    Returns:
        (ra_deg, dec_deg, size_arcsec) as floats, or None if pattern does not match.
    """
    base = os.path.basename(filename)
    name, _ = os.path.splitext(base)

    # Regex to capture RA, Dec, size
    pattern = r"_RA([0-9.+-]+)_Dec([0-9.+-]+)_size([0-9.+-]+)$"
    m = re.search(pattern, name)
    if m is None:
        return None

    ra_str, dec_str, size_str = m.groups()
    try:
        ra_deg = float(ra_str)
        dec_deg = float(dec_str)
        size_arcsec = float(size_str)
    except ValueError:
        return None

    return ra_deg, dec_deg, size_arcsec


def gaia_stars_for_cutout(
    ra_center,
    dec_center,
    size_px,
    size_arcsec,
    gaia_tab,
    gaia_coord,
    pixscale,
    maglim,
):
    """Return an Astropy Table of Gaia stars inside a square cutout.

    Columns included:
      ra_deg, dec_deg, Gmag, sep_arcsec, sep_pix, x_pix, y_pix, (optional) source_id
    """
    # Center coordinate
    center = SkyCoord(ra_center * u.deg, dec_center * u.deg)

    # Half diagonal of the square cutout in arcsec
    half_diag_arcsec = 0.5 * np.sqrt(2.0) * size_arcsec

    # Compute separation from center to all Gaia stars
    sep2d = center.separation(gaia_coord)

    # Find stars within the search radius
    idx_gaia = np.where(sep2d < half_diag_arcsec * u.arcsec)[0]

    if len(idx_gaia) == 0:
        return None

    # Magnitude cut
    mag = gaia_tab["phot_g_mean_mag"][idx_gaia]
    m = mag < maglim
    if not np.any(m):
        return None

    idx_gaia = idx_gaia[m]
    sep2d = sep2d[idx_gaia]
    mag = mag[m]

    # Small-angle conversion to pixel coordinates
    ra0 = ra_center
    dec0 = dec_center
    dec0_rad = np.deg2rad(dec0)

    dra_deg = gaia_tab["ra"][idx_gaia] - ra0
    ddec_deg = gaia_tab["dec"][idx_gaia] - dec0

    # Tangent-plane offsets in arcsec (east, north)
    dra_arcsec = dra_deg * np.cos(dec0_rad) * 3600.0   # +east on sky
    ddec_arcsec = ddec_deg * 3600.0                    # +north

    # Pixel offsets: x increasing to the right (west), y up (north)
    dx_pix = -dra_arcsec / pixscale
    dy_pix =  ddec_arcsec / pixscale

    x_pix = size_px / 2.0 + dx_pix
    y_pix = size_px / 2.0 + dy_pix

    # Keep only stars whose centers fall inside the cutout
    inside = (
        (x_pix >= 0.0) & (x_pix < size_px) &
        (y_pix >= 0.0) & (y_pix < size_px)
    )

    if not np.any(inside):
        return None

    x_pix = x_pix[inside]
    y_pix = y_pix[inside]
    mag = mag[inside]
    idx_gaia = idx_gaia[inside]
    sep_arcsec = sep2d[inside].arcsec
    sep_pix = sep_arcsec / pixscale

    # Build output table
    out = Table()
    if "source_id" in gaia_tab.colnames:
        out["source_id"] = gaia_tab["source_id"][idx_gaia]

    out["ra_deg"] = gaia_tab["ra"][idx_gaia]
    out["dec_deg"] = gaia_tab["dec"][idx_gaia]
    out["Gmag"] = mag
    out["sep_arcsec"] = sep_arcsec
    out["sep_pix"] = sep_pix
    out["x_pix"] = x_pix
    out["y_pix"] = y_pix

    return out


def write_gaia_sidecar(sidecar_path_base, gaia_table, overwrite=False):
    """Write Gaia star table for a cutout as a human-readable CSV sidecar."""
    if gaia_table is None or len(gaia_table) == 0:
        return

    out_path = sidecar_path_base + "_gaia.csv"
    if os.path.exists(out_path) and not overwrite:
        return

    gaia_table.write(out_path, format="ascii.csv", overwrite=True)
    return out_path


def main():
    args = parse_args()

    cutout_folder = args.cutout_folder
    pixscale = args.pixscale

    if not os.path.isdir(cutout_folder):
        raise NotADirectoryError(f"Cutout folder not found: {cutout_folder}")

    gaia_tab, gaia_coord = load_gaia_catalog(args.gaia_catalog)

    # List JPEG files in folder
    files = sorted(
        f for f in os.listdir(cutout_folder)
        if f.lower().endswith((".jpeg", ".jpg"))
    )

    print(f"Found {len(files)} JPEG cutouts in {cutout_folder}")
    print(f"Using Gaia magnitude limit: g < {args.gaia_maglim}")
    print()

    n_sidecars = 0
    n_with_stars = 0

    for i, fname in enumerate(files):
        full_path = os.path.join(cutout_folder, fname)
        parsed = parse_cutout_filename(fname)
        if parsed is None:
            continue

        ra_c, dec_c, size_arcsec = parsed
        size_px = int(round(size_arcsec / pixscale))

        gaia_tbl = gaia_stars_for_cutout(
            ra_center=ra_c,
            dec_center=dec_c,
            size_px=size_px,
            size_arcsec=size_arcsec,
            gaia_tab=gaia_tab,
            gaia_coord=gaia_coord,
            pixscale=pixscale,
            maglim=args.gaia_maglim,
        )

        if gaia_tbl is None or len(gaia_tbl) == 0:
            continue

        # Base path without extension, for sidecar naming
        base_name, _ = os.path.splitext(full_path)
        out_path = write_gaia_sidecar(
            sidecar_path_base=base_name,
            gaia_table=gaia_tbl,
            overwrite=args.overwrite,
        )

        if out_path:
            n_sidecars += 1
            n_with_stars += 1
            if (n_sidecars % 1000) == 0:
                print(f"Processed {i+1}/{len(files)} cutouts, wrote {n_sidecars} sidecars")

    print()
    print(f"Finished! Wrote {n_sidecars} sidecar files for {n_with_stars} cutouts with Gaia stars.")


if __name__ == "__main__":
    main()
