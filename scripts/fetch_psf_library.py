#!/usr/bin/env python
"""Fetch SGA-2020 PSF stamps and calibrate a parametric PSF profile.

SGA-2020 ships a pixelised grz PSF per galaxy group. These stamps are a
one-off calibration set, not a runtime library: they exist to fit the
parametric profile shape a forward model can use, and to measure how
well a kernel built from a single FWHM number reproduces a real coadd
PSF — the consumer reads that FWHM per sample from the PSF_G/R/Z header
cards, so this is the accuracy ceiling of its PSF model.

Run scripts/fetch_catalogues.py once first. For --n-groups groups
spread over the RA range of the selection:

  1. fetch the three per-band stamps,
  2. fit a free circular Moffat (amplitude, centre, alpha, beta) to
     each stamp,
  3. rebuild each stamp from its fitted FWHM alone, holding beta at the
     per-band median — exactly what a consumer of the PSF_* header cards
     can do — and report the total-variation residual; a Gaussian at the
     same FWHM is reported alongside to show what the profile wings are
     worth.

Fitted parameters land in $SGADATA/psf/psf_calibration.json, a summary
figure in $SGADATA/psf/psf_calibration.png.

Example:
    python scripts/fetch_psf_library.py --n-groups 300
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import fetch, urls
from src.bricks import Bricks
from src.catalog import Catalog
from src.cutout import BANDS, PIXSCALE
from src.select import galaxy_targets


def moffat(params, xs, ys):
    amplitude, x0, y0, alpha, beta = params
    r2 = (xs - x0) ** 2 + (ys - y0) ** 2
    return amplitude * (1.0 + r2 / alpha ** 2) ** (-beta)


def moffat_fwhm(alpha, beta):
    return 2.0 * alpha * np.sqrt(2.0 ** (1.0 / beta) - 1.0)


def fit_moffat(stamp):
    """Free circular-Moffat fit; returns (params, fwhm_px) or None."""
    ny, nx = stamp.shape
    ys, xs = np.mgrid[0:ny, 0:nx].astype(float)
    peak = float(stamp.max())
    if not peak > 0:
        return None
    y0, x0 = np.unravel_index(np.argmax(stamp), stamp.shape)
    start = [peak, float(x0), float(y0), 3.0, 2.5]
    bounds = ([0.0, x0 - 3.0, y0 - 3.0, 0.3, 1.1],
              [10.0 * peak, x0 + 3.0, y0 + 3.0, 30.0, 10.0])
    result = least_squares(
        lambda p: (moffat(p, xs, ys) - stamp).ravel(), start, bounds=bounds
    )
    if not result.success:
        return None
    return result.x, float(moffat_fwhm(result.x[3], result.x[4]))


def kernel_from_fwhm(shape, x0, y0, fwhm_px, beta=None):
    """Normalised kernel at a given FWHM: Moffat with fixed beta, or Gaussian."""
    ny, nx = shape
    ys, xs = np.mgrid[0:ny, 0:nx].astype(float)
    r2 = (xs - x0) ** 2 + (ys - y0) ** 2
    if beta is None:
        sigma = fwhm_px / 2.3548
        kernel = np.exp(-0.5 * r2 / sigma ** 2)
    else:
        alpha = fwhm_px / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))
        kernel = (1.0 + r2 / alpha ** 2) ** (-beta)
    return kernel / kernel.sum()


def total_variation(stamp, kernel):
    """Half the L1 distance between two unit-sum kernels, in [0, 1]."""
    return 0.5 * float(np.abs(stamp - kernel).sum())


def pick_groups(args):
    targets = galaxy_targets(Catalog(), Bricks(), zcut=args.zcut)
    groups = sorted(
        {(t["group_name"], t["group_ra"]) for t in targets},
        key=lambda g: g[1],
    )
    if args.n_groups < len(groups):
        picks = np.linspace(0, len(groups) - 1, args.n_groups).round().astype(int)
        groups = [groups[i] for i in np.unique(picks)]
    return groups


def main(args):
    psf_dir = Path(args.root) / "psf"
    psf_dir.mkdir(exist_ok=True, parents=True)
    groups = pick_groups(args)
    print(f"{len(groups)} groups")

    tasks = []
    for group_name, group_ra in groups:
        for name, url in urls.psf_stamps(group_ra, group_name).items():
            tasks.append((url, psf_dir / name))
    transferred, failures = fetch.fetch_many(tasks, workers=args.workers,
                                             desc="psf stamps")
    print(f"{transferred / 1e6:.0f} MB fetched, "
          f"{len(failures)} unavailable (groups without a largegalaxy run)")

    fits_by_band = {band: [] for band in BANDS}
    for group_name, _ in groups:
        for band in BANDS:
            path = psf_dir / f"{group_name}-largegalaxy-psf-{band}.fits.fz"
            if not path.exists():
                continue
            with fits.open(path) as hdul:
                hdu = hdul[1] if hdul[0].data is None else hdul[0]
                stamp = np.asarray(hdu.data, dtype=float)
            stamp /= stamp.sum()
            fitted = fit_moffat(stamp)
            if fitted is None:
                continue
            (amplitude, x0, y0, alpha, beta), fwhm = fitted
            fits_by_band[band].append({
                "group": group_name, "alpha_px": alpha, "beta": beta,
                "fwhm_px": fwhm, "x0": x0, "y0": y0,
                "stamp_path": str(path),
            })

    summary = {}
    residual_curves = {}
    for band in BANDS:
        rows = fits_by_band[band]
        betas = np.array([r["beta"] for r in rows])
        beta_band = float(np.median(betas))
        moffat_res, gauss_res = [], []
        for r in rows:
            with fits.open(r["stamp_path"]) as hdul:
                hdu = hdul[1] if hdul[0].data is None else hdul[0]
                stamp = np.asarray(hdu.data, dtype=float)
            stamp /= stamp.sum()
            moffat_res.append(total_variation(stamp, kernel_from_fwhm(
                stamp.shape, r["x0"], r["y0"], r["fwhm_px"], beta=beta_band)))
            gauss_res.append(total_variation(stamp, kernel_from_fwhm(
                stamp.shape, r["x0"], r["y0"], r["fwhm_px"])))
        moffat_res, gauss_res = np.array(moffat_res), np.array(gauss_res)
        residual_curves[band] = (moffat_res, gauss_res)
        summary[band] = {
            "n_stamps": len(rows),
            "beta_median": beta_band,
            "beta_quartiles": [float(np.percentile(betas, 25)),
                               float(np.percentile(betas, 75))],
            "fwhm_median_arcsec": float(np.median(
                [r["fwhm_px"] for r in rows]) * PIXSCALE),
            "tv_moffat_median": float(np.median(moffat_res)),
            "tv_moffat_p90": float(np.percentile(moffat_res, 90)),
            "tv_gauss_median": float(np.median(gauss_res)),
            "tv_gauss_p90": float(np.percentile(gauss_res, 90)),
        }
        print(f"{band}: {len(rows)} stamps | beta median {beta_band:.2f} "
              f"[{summary[band]['beta_quartiles'][0]:.2f}, "
              f"{summary[band]['beta_quartiles'][1]:.2f}] | "
              f"FWHM median {summary[band]['fwhm_median_arcsec']:.2f}\" | "
              f"TV residual: Moffat {summary[band]['tv_moffat_median']:.3f} "
              f"(p90 {summary[band]['tv_moffat_p90']:.3f}), "
              f"Gaussian {summary[band]['tv_gauss_median']:.3f}")

    with open(psf_dir / "psf_calibration.json", "w") as f:
        json.dump({"summary": summary, "fits": fits_by_band}, f, indent=1)
    print(f"wrote {psf_dir / 'psf_calibration.json'}")
    figure(psf_dir, fits_by_band, residual_curves, summary)


def figure(psf_dir, fits_by_band, residual_curves, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for k, band in enumerate(BANDS):
        rows = fits_by_band[band]
        beta_band = summary[band]["beta_median"]
        axes[0, k].hist([r["beta"] for r in rows], bins=30, color="C0")
        axes[0, k].axvline(beta_band, color="C3")
        axes[0, k].set_title(f"{band}: Moffat beta (median {beta_band:.2f})")
        moffat_res, gauss_res = residual_curves[band]
        bins = np.linspace(0, max(gauss_res.max(), 0.1), 40)
        axes[1, k].hist(moffat_res, bins=bins, alpha=0.7,
                        label="Moffat @ median beta")
        axes[1, k].hist(gauss_res, bins=bins, alpha=0.7,
                        label="Gaussian")
        axes[1, k].set_title(f"{band}: TV residual of FWHM-only kernel")
        axes[1, k].legend()
    fig.tight_layout()
    out = psf_dir / "psf_calibration.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--root", default="data/sga",
                        help="data root: PSF stamps land under <root>/psf")
    parser.add_argument("--zcut", type=float, default=0.05,
                        help="upper redshift cut on Z_LEDA")
    parser.add_argument("--n-groups", type=int, default=300,
                        help="groups to fetch, spread over the RA range")
    parser.add_argument("--workers", type=int, default=12,
                        help="parallel download connections")
    main(parser.parse_args())
