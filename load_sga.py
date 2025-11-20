#! /usr/bin/env python3

import os
import sys
import argparse
import math

import numpy as np
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from tqdm.contrib.concurrent import process_map


parser = argparse.ArgumentParser(description="Download SGA galaxy images")
parser.add_argument('--size', type=int, default='-1', help='Force image size, -1 is an adaptive resolution based on the actual galaxy size')
parser.add_argument('--catalog', type=str, default='./data/catalogues/SGA-2020.fits', help='Path to the SGA catalogue file')
parser.add_argument('--zcut', type=float, default=0.05, help='Redshift cutoff for selecting galaxies')
parser.add_argument('--na', type=float, default=3.0, help='Number of semi-major axes to cover if using adaptive size')
parser.add_argument('--pixscale', type=float, default=0.262, help='Pixel scale in arcseconds per pixel')
parser.add_argument('--layer', type=str, default='ls-dr9', help='Layer to use for the image cutout')
parser.add_argument('--download_folder', type=str, default='./data/SGA', help='Folder to save downloaded images')
parser.add_argument('--num_workers', type=int, default=4, help='Number of parallel workers for downloading images')
parser.add_argument('--chunksize', type=int, default=32, help='Chunk size for parallel processing')
parser.add_argument('--n_limit', type=int, default=None, help='Limit to N galaxies (for testing)')
parser.add_argument('--mode', type=str, default='galaxy', choices=['galaxy', 'fields', 'retry'],
                    help='Download mode: galaxy (central), fields (background), or retry (redownload missing/broken cutouts)')
parser.add_argument('--n_away', type=int, default=7, help='Number of field positions per galaxy (mode=fields)')
parser.add_argument('--na_away', type=float, default=20.0, help='Offset multiplier for field radius (na_away * D26/2)')
parser.add_argument('--bricks_file', type=str, default='./data/catalogues/survey-bricks-dr9-south.fits',
                    help='Path to survey-bricks-dr9-south.fits[.gz] for footprint checking (default: ./data/catalogues/survey-bricks-dr9-south.fits)')
parser.add_argument('--bricks_file_north', type=str, default='./data/catalogues/survey-bricks-dr9-north.fits',
                    help='Path to survey-bricks-dr9-north.fits[.gz] for northern footprint (default: ./data/catalogues/survey-bricks-dr9-north.fits)')
parser.add_argument('--min_brick_exposures', type=int, default=1,
                    help='Minimum total exposures (g+r+z) in brick to consider it covered')
parser.add_argument('--broken_size_threshold', type=int, default=500,
                    help='File size threshold (bytes) to identify broken downloads (HTML errors)')
args = parser.parse_args()


def load_decals_bricks(bricks_path):
    """Load DECaLS DR9 bricks and build SkyCoord of centers."""
    if bricks_path is None:
        return None, None

    if not os.path.exists(bricks_path):
        raise FileNotFoundError(f"Bricks file not found: {bricks_path}")

    with fits.open(bricks_path) as hdul:
        bricks = hdul[1].data
        brick_centers = SkyCoord(bricks['ra'] * u.deg, bricks['dec'] * u.deg)

    print(f"Loaded {len(bricks)} DECaLS DR9 bricks from {os.path.basename(bricks_path)}")
    return bricks, brick_centers


def in_decals_footprint(ra_deg, dec_deg, bricks_list, brick_centers_list, min_exposures=1):
    """Vectorized check: is each galaxy in DECaLS DR9 footprint (north or south)?

    Parameters:
    - bricks_list: list of bricks arrays (can contain south, north, or both)
    - brick_centers_list: list of SkyCoord objects for brick centers
    - min_exposures: minimum total exposures (g+r+z) required
    """
    ra_deg = np.asarray(ra_deg)
    dec_deg = np.asarray(dec_deg)

    coords = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)

    # Initialize mask as False for all galaxies
    in_footprint = np.zeros(len(ra_deg), dtype=bool)

    # Check against each bricks catalog (south and/or north)
    for bricks, brick_centers in zip(bricks_list, brick_centers_list):
        # Find nearest brick for each galaxy
        idx_brick, sep2d, _ = coords.match_to_catalog_sky(brick_centers)

        # Check exposure count
        if min_exposures is not None:
            nexp = (bricks['nexp_g'][idx_brick] +
                    bricks['nexp_r'][idx_brick] +
                    bricks['nexp_z'][idx_brick])
            m_exp = nexp >= min_exposures
        else:
            m_exp = np.ones_like(idx_brick, dtype=bool)

        # Get brick bounds
        ra1 = bricks['ra1'][idx_brick]
        ra2 = bricks['ra2'][idx_brick]
        dec1 = bricks['dec1'][idx_brick]
        dec2 = bricks['dec2'][idx_brick]

        # RA check with wrap-around
        width = (ra2 - ra1 + 360.0) % 360.0
        dra = (ra_deg - ra1 + 360.0) % 360.0
        m_ra = (dra >= 0.0) & (dra <= width)

        # Dec check
        m_dec = (dec_deg >= dec1) & (dec_deg <= dec2)

        # Mark galaxies in this catalog as in footprint
        in_footprint |= (m_exp & m_ra & m_dec)

    return in_footprint


def generate_field_positions(ra_center, dec_center, d26, n_away, na_away, pixscale):
    """Generate uniformly distributed field positions around a galaxy.

    Parameters:
    - ra_center, dec_center: galaxy center in degrees
    - d26: diameter at 26 mag/arcsec² in arcmin
    - n_away: number of field positions
    - na_away: offset multiplier (radius = na_away × D26/2 in arcsec)
    - pixscale: pixel scale in arcsec/pixel

    Returns:
    - List of (ra, dec, angle_deg, size_px, size_arcsec) tuples
    """
    # Radius in arcsec: na_away × (D26/2 arcmin)
    # Factor 30 = 60 (arcmin→arcsec) / 2 (diameter→radius)
    radius_arcsec = na_away * d26 * 30.0
    radius_deg = radius_arcsec / 3600.0

    # Size in pixels: na=3.0 (default) × (D26/2 in arcsec) / pixscale
    size_px = int(max(16, min(512, round(3.0 * d26 * 30.0 / pixscale))))
    size_arcsec = 3.0 * d26 * 30.0

    positions = []
    for i in range(n_away):
        # Angle in radians and degrees
        angle_rad = 2 * math.pi * i / n_away
        angle_deg = math.degrees(angle_rad)

        # Account for declination when computing RA offset
        dec_rad = math.radians(dec_center)
        ra_offset = radius_deg / math.cos(dec_rad) if abs(dec_rad) < math.pi / 2 else radius_deg
        dec_offset = radius_deg

        # Position on the circle
        ra = ra_center + ra_offset * math.cos(angle_rad)
        dec = dec_center + dec_offset * math.sin(angle_rad)

        # Clamp Dec to valid range
        dec = max(-90, min(90, dec))

        positions.append((ra, dec, angle_deg, size_px, size_arcsec))

    return positions


def download_single(ra, dec, idx, size_px, size_arcsec, angle=None):
    """Download a single cutout (galaxy or field).

    Parameters:
    - ra, dec: coordinates in degrees
    - idx: SGA galaxy ID
    - size_px: cutout size in pixels
    - size_arcsec: cutout size in arcseconds
    - angle: position angle for field cutout (None for central galaxy)
    """
    # Prepare the parameters for the request
    params = {'ra': ra, 'dec': dec, 'size': size_px, 'layer': args.layer, 'pixscale': args.pixscale}

    # Create the full URL
    url = requests.Request('GET', 'https://www.legacysurvey.org/viewer/jpeg-cutout', params=params).prepare().url

    # Request and download the image
    response = requests.get(url)

    # Define filename based on mode
    if angle is None:
        # Central galaxy
        filename = f"sga_{idx}_RA{ra:.4f}_Dec{dec:.4f}_size{size_arcsec:.1f}.jpeg"
    else:
        # Field cutout with position angle
        filename = f"sga_{idx}_field_{angle:.1f}_RA{ra:.4f}_Dec{dec:.4f}_size{size_arcsec:.1f}.jpeg"

    file_path = os.path.join(args.download_folder, filename)

    # Save the downloaded content to a file
    with open(file_path, "wb") as f:
        f.write(response.content)


def download_wrapper(download_tuple):
    """Wrapper to unpack tuple and call download_single (for multiprocessing)."""
    return download_single(*download_tuple)


def find_missing_cutouts(download_folder, catalog_path, zcut, broken_size_threshold, n_limit=None):
    """Identify broken or missing cutouts in a folder.

    Parameters:
    - download_folder: folder containing downloaded JPEGs
    - catalog_path: path to SGA catalog
    - zcut: redshift cutoff
    - broken_size_threshold: max file size to consider broken (HTML error)
    - n_limit: limit to first N galaxies (optional)

    Returns:
    - List of (ra, dec, sga_id, d26) tuples for missing/broken files
    """
    # Load catalog
    data = fits.open(catalog_path)[1].data
    sel = (data['Z_LEDA'] > 0) & (data['Z_LEDA'] < zcut)

    ra_arr = data['RA'][sel]
    dec_arr = data['DEC'][sel]
    id_arr = data['SGA_ID'][sel]
    d26_arr = data['D26'][sel]

    # Apply limit if specified
    if n_limit is not None:
        ra_arr = ra_arr[:n_limit]
        dec_arr = dec_arr[:n_limit]
        id_arr = id_arr[:n_limit]
        d26_arr = d26_arr[:n_limit]

    # Find all JPEG files in download folder
    if not os.path.exists(download_folder):
        print(f"Download folder {download_folder} not found")
        return []

    existing_files = {}
    for filename in os.listdir(download_folder):
        if filename.endswith('.jpeg'):
            # Parse SGA_ID from filename (format: sga_XXXX_...)
            try:
                sga_id = int(filename.split('_')[1])
                filepath = os.path.join(download_folder, filename)
                file_size = os.path.getsize(filepath)
                existing_files[sga_id] = (filepath, file_size)
            except (ValueError, IndexError):
                pass

    # Identify broken/missing
    missing = []
    broken = []

    for i in range(len(id_arr)):
        sga_id = int(id_arr[i])
        if sga_id not in existing_files:
            # Missing file
            missing.append((ra_arr[i], dec_arr[i], sga_id, d26_arr[i]))
        else:
            # Check if file is broken (HTML error, small file size)
            filepath, file_size = existing_files[sga_id]
            if file_size < broken_size_threshold:
                broken.append((ra_arr[i], dec_arr[i], sga_id, d26_arr[i]))

    print(f"Status in {download_folder}:")
    print(f"  Total galaxies: {len(id_arr)}")
    print(f"  Downloaded: {len(existing_files)}")
    print(f"  Missing files: {len(missing)}")
    print(f"  Broken files: {len(broken)}")
    print(f"  Total to retry: {len(missing) + len(broken)}\n")

    return missing + broken


if __name__ == '__main__':
    # Handle retry mode: identify and redownload missing/broken cutouts
    if args.mode == 'retry':
        print(f"Retry mode: checking {args.download_folder}\n")
        missing_cutouts = find_missing_cutouts(
            args.download_folder,
            args.catalog,
            args.zcut,
            args.broken_size_threshold,
            n_limit=args.n_limit
        )

        if not missing_cutouts:
            print("No missing or broken files found!")
            sys.exit(0)

        # Prepare download list for missing cutouts
        download_list = []
        for ra, dec, sga_id, d26 in missing_cutouts:
            if args.size == -1:
                size_px = int(max(16, min(512, round(args.na * d26 * 30.0 / args.pixscale))))
                size_arcsec = args.na * d26 * 30.0
            else:
                size_px = args.size
                size_arcsec = args.size * args.pixscale
            download_list.append((ra, dec, sga_id, size_px, size_arcsec, None))

        print(f"Attempting to redownload {len(download_list)} cutouts...\n")
        process_map(
            download_wrapper,
            download_list,
            total=len(download_list),
            desc='Retrying downloads',
            max_workers=args.num_workers,
            chunksize=args.chunksize
        )
        print(f"Retry complete!")
        sys.exit(0)

    # Normal mode: load catalog and download
    os.makedirs(args.download_folder, exist_ok=True)
    data = fits.open(args.catalog)[1].data

    sel = (data['Z_LEDA'] > 0) & (data['Z_LEDA'] < args.zcut)
    nsel = sel.sum()
    print(f"Number of objects with Z_LEDA < {args.zcut}: {nsel}")


    ra_arr = data['RA'][sel]
    dec_arr = data['DEC'][sel]
    id_arr = data['SGA_ID'][sel]
    d26_arr = data['D26'][sel]

    # Apply limit if specified
    if args.n_limit is not None:
        ra_arr = ra_arr[:args.n_limit]
        dec_arr = dec_arr[:args.n_limit]
        id_arr = id_arr[:args.n_limit]
        d26_arr = d26_arr[:args.n_limit]
        nsel = min(nsel, args.n_limit)

    # Optional: restrict to DECaLS footprint using bricks (north and/or south)
    if args.bricks_file is not None or args.bricks_file_north is not None:
        bricks_list = []
        brick_centers_list = []

        # Load southern bricks if provided
        if args.bricks_file is not None:
            bricks_s, brick_centers_s = load_decals_bricks(args.bricks_file)
            if bricks_s is None:
                raise RuntimeError("Could not load bricks file, but --bricks_file was given.")
            bricks_list.append(bricks_s)
            brick_centers_list.append(brick_centers_s)

        # Load northern bricks if provided
        if args.bricks_file_north is not None:
            bricks_n, brick_centers_n = load_decals_bricks(args.bricks_file_north)
            if bricks_n is None:
                raise RuntimeError("Could not load northern bricks file, but --bricks_file_north was given.")
            bricks_list.append(bricks_n)
            brick_centers_list.append(brick_centers_n)

        in_fp = in_decals_footprint(
            ra_arr,
            dec_arr,
            bricks_list,
            brick_centers_list,
            min_exposures=args.min_brick_exposures,
        )

        n_in = in_fp.sum()
        n_out = nsel - n_in
        region_str = "DECaLS DR9 (north+south)" if (args.bricks_file and args.bricks_file_north) else ("DECaLS DR9-north" if args.bricks_file_north else "DECaLS DR9-south")
        print(f"{n_in} of {nsel} galaxies lie in {region_str} footprint")
        if n_out > 0:
            print(f"  ({n_out} galaxies filtered out)")

        # Apply mask to all arrays
        ra_arr = ra_arr[in_fp]
        dec_arr = dec_arr[in_fp]
        id_arr = id_arr[in_fp]
        d26_arr = d26_arr[in_fp]
        nsel = n_in

    if args.size == -1:
        # Adaptive sizing: size = na × (D26/2 in arcsec) / pixscale
        # D26 is in arcmin (diameter), so: * 30 = * 60 (arcmin→arcsec) / 2 (diameter→radius)
        size_px = (args.na * d26_arr * 30.0 / args.pixscale).astype(int)
        size_arcsec = args.na * d26_arr * 30.0
    else:
        size_px = np.full(nsel, args.size)
        size_arcsec = np.full(nsel, args.size * args.pixscale)

    # Prepare download list based on mode
    if args.mode == 'galaxy':
        # Download central galaxies only
        download_list = []
        for i in range(nsel):
            download_list.append((ra_arr[i], dec_arr[i], id_arr[i], size_px[i], size_arcsec[i], None))

        process_map(
            download_wrapper,
            download_list,
            total=len(download_list),
            desc=f'Downloading galaxy images',
            max_workers=args.num_workers,
            chunksize=args.chunksize
        )
        print(f"Downloaded {nsel} galaxy images to {args.download_folder}")

    elif args.mode == 'fields':
        # Download field cutouts only (no central galaxies)
        download_list = []
        for i in range(nsel):
            field_positions = generate_field_positions(
                ra_arr[i], dec_arr[i], d26_arr[i],
                args.n_away, args.na_away, args.pixscale
            )
            for ra, dec, angle, size_px_field, size_arcsec_field in field_positions:
                download_list.append((ra, dec, id_arr[i], size_px_field, size_arcsec_field, angle))

        n_total = len(download_list)
        process_map(
            download_wrapper,
            download_list,
            total=n_total,
            desc=f'Downloading field images',
            max_workers=args.num_workers,
            chunksize=args.chunksize
        )
        print(f"Downloaded {n_total} field images ({args.n_away} per galaxy × {nsel} galaxies) to {args.download_folder}")
