#! /usr/bin/env python3

import os
import re
import argparse

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create Gaia bright-star sidecar CSV files for DECALS cutouts (global search)."
    )
    parser.add_argument(
        "--cutout_folder",
        type=str,
        required=True,
        help="Folder containing JPEG cutouts (from your SGA downloader).",
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
        help="Pixel scale in arcsec per pixel (must match what you used for downloading).",
    )
    parser.add_argument(
        "--gaia_maglim",
        type=float,
        default=17.0,
        help="G-band magnitude limit for Gaia bright stars.",
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
    Optionally: source_id.

    Loads only necessary columns to save memory (filters out parallax, proper motion, etc.).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Gaia catalog file not found: {path}")

    # Read full catalog first to check structure
    gaia = Table.read(path)
    required_cols = ["ra", "dec", "phot_g_mean_mag"]
    if not all(col in gaia.colnames for col in required_cols):
        raise ValueError(
            "Gaia catalog must contain columns: ra, dec, phot_g_mean_mag"
        )

    # Select only essential columns for this task
    cols_to_keep = required_cols.copy()
    if "source_id" in gaia.colnames:
        cols_to_keep.insert(0, "source_id")

    gaia = gaia[cols_to_keep]

    gaia_coord = SkyCoord(gaia["ra"] * u.deg, gaia["dec"] * u.deg)
    print(f"Loaded Gaia catalog with {len(gaia)} entries ({len(gaia.colnames)} columns) from {path}")
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


def write_gaia_sidecar(sidecar_path_base, gaia_table, overwrite=False):
    """Write Gaia star table for a cutout as a human-readable CSV sidecar.

    sidecar_path_base should be the full path without extension.
    """
    if gaia_table is None or len(gaia_table) == 0:
        return

    out_path = sidecar_path_base + "_gaia.csv"
    if os.path.exists(out_path) and not overwrite:
        return

    gaia_table.write(out_path, format="ascii.csv", overwrite=True)


def main():
    args = parse_args()

    cutout_folder = args.cutout_folder
    pixscale = args.pixscale

    if not os.path.isdir(cutout_folder):
        raise NotADirectoryError(f"Cutout folder not found: {cutout_folder}")

    # Load Gaia catalog and coordinates
    gaia_tab, gaia_coord = load_gaia_catalog(args.gaia_catalog)

    # Collect cutouts and parse their metadata
    jpeg_files = sorted(
        f for f in os.listdir(cutout_folder)
        if f.lower().endswith((".jpeg", ".jpg"))
    )

    cutout_paths = []
    base_paths = []
    ra_centers = []
    dec_centers = []
    size_arcsec_list = []

    for fname in jpeg_files:
        full_path = os.path.join(cutout_folder, fname)
        parsed = parse_cutout_filename(fname)
        if parsed is None:
            continue

        ra_c, dec_c, size_arcsec = parsed

        cutout_paths.append(full_path)
        base_paths.append(os.path.splitext(full_path)[0])
        ra_centers.append(ra_c)
        dec_centers.append(dec_c)
        size_arcsec_list.append(size_arcsec)

    if len(cutout_paths) == 0:
        print("No parsable cutout files found. Nothing to do.")
        return

    ra_centers = np.array(ra_centers, dtype=float)
    dec_centers = np.array(dec_centers, dtype=float)
    size_arcsec_arr = np.array(size_arcsec_list, dtype=float)

    # Approximate size in pixels from size_arcsec and pixscale
    # (may differ by 1 pixel from downloader int(), but fine for star positions)
    size_px_arr = np.rint(size_arcsec_arr / pixscale).astype(int)

    n_cutouts = len(cutout_paths)
    print(f"Found {n_cutouts} cutouts with parsable metadata.")

    # Build SkyCoord for all cutout centers
    centers = SkyCoord(ra_centers * u.deg, dec_centers * u.deg)

    # Per cutout half-diagonal in arcsec and global max radius
    half_diag_arcsec = 0.5 * np.sqrt(2.0) * size_arcsec_arr
    max_radius = half_diag_arcsec.max()
    print(f"Global search radius (max half diagonal) = {max_radius:.2f} arcsec")

    # Single global search: all Gaia stars around all cutout centers
    # gaia_coord.search_around_sky(centers, ...) returns:
    #   idx_center: indices into centers (cutouts)
    #   idx_gaia:   indices into gaia_coord (Gaia rows)
    idx_center, idx_gaia, sep2d, _ = gaia_coord.search_around_sky(
        centers, max_radius * u.arcsec
    )

    print(f"Global search found {len(idx_gaia)} (cutout, Gaia) candidate pairs.")

    if len(idx_gaia) == 0:
        print("No Gaia stars found within max_radius of any cutout. Done.")
        return

    # Sort matches by cutout index so we can slice per cutout efficiently
    order = np.argsort(idx_center)
    idx_center_s = idx_center[order]
    idx_gaia_s = idx_gaia[order]
    sep2d_s = sep2d[order]

    unique_centers, start_indices = np.unique(idx_center_s, return_index=True)
    # add end index sentinel
    end_indices = np.empty_like(start_indices)
    end_indices[:-1] = start_indices[1:]
    end_indices[-1] = len(idx_center_s)

    n_sidecars = 0

    # Loop over only the cutouts that actually have at least one candidate match
    for k, c_idx in enumerate(unique_centers):
        i = int(c_idx)  # cutout index

        row_start = int(start_indices[k])
        row_end = int(end_indices[k])

        gi_block = idx_gaia_s[row_start:row_end]
        sep_block_arcsec = sep2d_s[row_start:row_end].arcsec

        # Enforce cutout specific radius (half diagonal) inside the global max
        r_i = half_diag_arcsec[i]
        m_radius = sep_block_arcsec <= r_i
        if not np.any(m_radius):
            continue

        gi = gi_block[m_radius]
        sep_arcsec = sep_block_arcsec[m_radius]

        # Magnitude cut
        mag = gaia_tab["phot_g_mean_mag"][gi]
        m_bright = mag < args.gaia_maglim
        if not np.any(m_bright):
            continue

        gi = gi[m_bright]
        sep_arcsec = sep_arcsec[m_bright]
        mag = mag[m_bright]

        # Small-angle conversion to pixel coordinates for this cutout
        ra0 = ra_centers[i]
        dec0 = dec_centers[i]
        size_px = size_px_arr[i]
        size_arcsec = size_arcsec_arr[i]

        dec0_rad = np.deg2rad(dec0)

        dra_deg = gaia_tab["ra"][gi] - ra0
        ddec_deg = gaia_tab["dec"][gi] - dec0

        # Tangent-plane offsets in arcsec (east, north)
        dra_arcsec = dra_deg * np.cos(dec0_rad) * 3600.0   # +east on sky
        ddec_arcsec = ddec_deg * 3600.0                    # +north

        # Pixel offsets: x increasing to the right (west), y up (north)
        dx_pix = -dra_arcsec / pixscale
        dy_pix =  ddec_arcsec / pixscale

        x_pix = size_px / 2.0 + dx_pix
        y_pix = size_px / 2.0 + dy_pix

        # Keep only stars whose centers fall inside the cutout frame
        inside = (
            (x_pix >= 0.0) & (x_pix < size_px) &
            (y_pix >= 0.0) & (y_pix < size_px)
        )

        if not np.any(inside):
            continue

        x_pix = x_pix[inside]
        y_pix = y_pix[inside]
        mag_i = mag[inside]
        gi_inside = gi[inside]
        sep_arcsec_inside = sep_arcsec[inside]
        sep_pix_inside = sep_arcsec_inside / pixscale

        # Build output table for this cutout
        tbl = Table()
        if "source_id" in gaia_tab.colnames:
            tbl["source_id"] = gaia_tab["source_id"][gi_inside]

        tbl["ra_deg"] = gaia_tab["ra"][gi_inside]
        tbl["dec_deg"] = gaia_tab["dec"][gi_inside]
        tbl["Gmag"] = mag_i
        tbl["sep_arcsec"] = sep_arcsec_inside
        tbl["sep_pix"] = sep_pix_inside
        tbl["x_pix"] = x_pix
        tbl["y_pix"] = y_pix

        # Write sidecar CSV for this cutout
        sidecar_base = base_paths[i]
        write_gaia_sidecar(sidecar_base, tbl, overwrite=args.overwrite)
        n_sidecars += 1

    print(f"Wrote {n_sidecars} sidecar files.")
    print("Gaia sidecar generation finished.")


if __name__ == "__main__":
    main()
