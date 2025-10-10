#! /usr/bin/env python3

import os
import argparse

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
args = parser.parse_args()


def download_single(ra, dec, idx, size):

    # Prepare the parameters for the request
    params = {'ra': ra, 'dec': dec, 'size': size, 'layer': args.layer, 'pixscale': args.pixscale}

    # Create the full URL
    url = requests.Request('GET', 'https://www.legacysurvey.org/viewer/jpeg-cutout', params=params).prepare().url

    # Request and download the image
    response = requests.get(url)

    # Define a unique filename for each galaxy image
    filename = f"sga_{idx}_RA{ra:.4f}_Dec{dec:.4f}.jpeg"
    file_path = os.path.join(args.download_folder, filename)

    # Save the downloaded content to a file
    with open(file_path, "wb") as f:
        f.write(response.content)


if __name__ == '__main__':
    os.makedirs(args.download_folder, exist_ok=True)
    data = fits.open(args.catalog)[1].data

    sel = (data['Z_LEDA'] > 0) & (data['Z_LEDA'] < args.zcut)
    nsel = sel.sum()
    print(f"Number of objects with Z_LEDA < {args.zcut}: {nsel}")


    ra_arr = data['RA'][sel]
    dec_arr = data['DEC'][sel]
    id_arr = data['SGA_ID'][sel]

    if args.size == -1:
        size = (args.na * data['D26'][sel] * 30.0 / args.pixscale).astype(int)
    else:
        size = np.full(nsel, args.size)

    process_map(
        download_single,
        ra_arr, dec_arr, id_arr, size,
        total=nsel,
        desc=f'Downloading images',
        max_workers=args.num_workers,
        chunksize=args.chunksize
    )

    print(f"Downloaded {nsel} images to {args.download_folder}")
