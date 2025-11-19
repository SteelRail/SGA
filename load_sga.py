#! /usr/bin/env python3

import os
import argparse
import math

import numpy as np
import requests
from astropy.io import fits
from tqdm.contrib.concurrent import process_map


parser = argparse.ArgumentParser(description="Download SGA galaxy images")
parser.add_argument('--size', type=int, default='-1', help='Force image size, -1 is an adaptive resolution based on the actual galaxy size')
parser.add_argument('--catalog', type=str, default='./data/catalogues/SGA-2020.fits', help='Path to the SGA catalogue file')
parser.add_argument('--zcut', type=float, default=0.05, help='Redshift cutoff for selecting galaxies')
parser.add_argument('--na', type=float, default=3.0, help='Number of semi-major axes to cover if using adaptive size')
parser.add_argument('--pixscale', type=float, default=0.262, help='Pixel scale in arcseconds per pixel')
parser.add_argument('--layer', type=str, default='ls-dr9', help='Layer to use for the image cutout')
parser.add_argument('--download_folder', type=str, default='./data/SGA', help='Folder to save downloaded images')
parser.add_argument('--num_workers', type=int, default=16, help='Number of parallel workers for downloading images')
parser.add_argument('--chunksize', type=int, default=32, help='Chunk size for parallel processing')
parser.add_argument('--n_limit', type=int, default=None, help='Limit to N galaxies (for testing)')
parser.add_argument('--mode', type=str, default='galaxy', choices=['galaxy', 'fields'],
                    help='Download mode: galaxy (central), or fields (background around galaxies)')
parser.add_argument('--n_away', type=int, default=7, help='Number of field positions per galaxy (mode=fields)')
parser.add_argument('--na_away', type=float, default=20.0, help='Offset multiplier for field radius (na_away * D26/2)')
args = parser.parse_args()


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


if __name__ == '__main__':
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
