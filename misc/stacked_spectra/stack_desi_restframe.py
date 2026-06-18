#!/usr/bin/env python
"""Shift DESI spectra to rest frame and stack them with a simple mean.

Example
-------
python stack_desi_restframe.py \
    --spectra spectra-*.fits \
    --redshifts zbest.fits \
    --output stacked-restframe.fits

The redshift table must contain TARGETID and Z columns.  Each input spectra file
is read directly with fitsio, matched to the redshift table by TARGETID, shifted
to rest wavelength via lambda_rest = lambda_obs / (1 + z), interpolated onto a
common grid, and combined with an unweighted mean.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import fitsio
import numpy as np
from astropy.io import fits
from astropy.table import Table


DEFAULT_REST_WAVE_MIN = 1200.0
DEFAULT_REST_WAVE_MAX = 10000.0
DEFAULT_REST_WAVE_STEP = 0.1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stack DESI spectra after shifting them into the rest frame."
    )
    parser.add_argument(
        "--spectra",
        nargs="+",
        required=True,
        help="Input DESI spectra FITS files. Shell globs are okay.",
    )
    parser.add_argument(
        "--redshifts",
        required=True,
        help="Redshift table with TARGETID and Z columns, e.g. zbest or redrock.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output FITS file for the stacked rest-frame spectrum.",
    )
    parser.add_argument(
        "--targetids",
        nargs="+",
        type=np.int64,
        default=None,
        help="Optional TARGETID subset to include.",
    )
    parser.add_argument(
        "--rest-wave-min",
        type=float,
        default=DEFAULT_REST_WAVE_MIN,
        help=f"Minimum rest-frame wavelength in Angstrom. Default: {DEFAULT_REST_WAVE_MIN}",
    )
    parser.add_argument(
        "--rest-wave-max",
        type=float,
        default=DEFAULT_REST_WAVE_MAX,
        help=f"Maximum rest-frame wavelength in Angstrom. Default: {DEFAULT_REST_WAVE_MAX}",
    )
    parser.add_argument(
        "--rest-wave-step",
        type=float,
        default=DEFAULT_REST_WAVE_STEP,
        help=f"Rest-frame wavelength grid spacing in Angstrom. Default: {DEFAULT_REST_WAVE_STEP}",
    )
    parser.add_argument(
        "--cameras",
        nargs="+",
        default=None,
        help="Optional camera list to use, e.g. b r z. Default: all cameras in each file.",
    )
    parser.add_argument(
        "--min-contrib",
        type=int,
        default=1,
        help="Set stacked flux to NaN where fewer than this many pixels contribute.",
    )
    parser.add_argument(
        "--clobber",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )
    return parser.parse_args()


def read_redshifts(path: str | Path, targetids: list[int] | None) -> dict[int, float]:
    table = Table.read(path)
    names = {name.upper(): name for name in table.colnames}
    if "TARGETID" not in names or "Z" not in names:
        raise ValueError("redshift table must contain TARGETID and Z columns")

    targetid_col = names["TARGETID"]
    z_col = names["Z"]
    keep = np.isfinite(table[z_col]) & (table[z_col] > -1.0)

    if targetids is not None:
        requested = np.asarray(targetids, dtype=np.int64)
        keep &= np.isin(np.asarray(table[targetid_col], dtype=np.int64), requested)

    return {
        int(targetid): float(z)
        for targetid, z in zip(table[targetid_col][keep], table[z_col][keep])
    }


def make_wave_grid(wmin: float, wmax: float, dw: float) -> np.ndarray:
    if dw <= 0:
        raise ValueError("--rest-wave-step must be positive")
    if wmax <= wmin:
        raise ValueError("--rest-wave-max must be greater than --rest-wave-min")
    nbin = int(np.floor((wmax - wmin) / dw)) + 1
    return wmin + dw * np.arange(nbin, dtype=float)


def valid_flux_mask(
    flux: np.ndarray, ivar: np.ndarray | None, mask: np.ndarray | None
) -> np.ndarray:
    valid = np.isfinite(flux)
    if ivar is not None:
        valid &= np.isfinite(ivar) & (ivar > 0)
    if mask is not None:
        valid &= mask == 0
    return valid


def add_camera_to_spectrum_grid(
    rest_wave_grid: np.ndarray,
    obs_wave: np.ndarray,
    flux: np.ndarray,
    redshift: float,
    spectrum_flux_sum: np.ndarray,
    spectrum_flux_count: np.ndarray,
    ivar: np.ndarray | None = None,
    mask: np.ndarray | None = None,
) -> None:
    keep = valid_flux_mask(flux, ivar, mask)
    if np.count_nonzero(keep) < 2:
        return

    rest_wave = obs_wave[keep] / (1.0 + redshift)
    rest_flux = flux[keep]
    order = np.argsort(rest_wave)
    rest_wave = rest_wave[order]
    rest_flux = rest_flux[order]

    interp_flux = np.interp(
        rest_wave_grid, rest_wave, rest_flux, left=np.nan, right=np.nan
    )
    good = np.isfinite(interp_flux)
    spectrum_flux_sum[good] += interp_flux[good]
    spectrum_flux_count[good] += 1


def _fits_extnames(hdus: fitsio.FITS) -> set[str]:
    return {hdus[i].get_extname().upper() for i in range(1, len(hdus))}


def _read_image_rows(hdu, rows: np.ndarray) -> np.ndarray:
    try:
        return hdu[rows, :]
    except Exception:
        return hdu.read()[rows]


def read_desi_spectra_rows(
    path: str | Path,
    wanted_targetids: set[int],
    cameras: list[str] | None = None,
) -> tuple[np.ndarray, dict[str, dict[str, np.ndarray | None]]]:
    camera_names = None if cameras is None else {camera.lower() for camera in cameras}

    with fitsio.FITS(path) as hdus:
        extnames = _fits_extnames(hdus)
        if "FIBERMAP" not in extnames:
            raise ValueError(f"{path} does not contain a FIBERMAP HDU")

        fibermap = hdus["FIBERMAP"].read(columns=["TARGETID"])
        targetids = np.asarray(fibermap["TARGETID"], dtype=np.int64)
        rows = np.flatnonzero(np.isin(targetids, list(wanted_targetids)))
        if rows.size == 0:
            return np.array([], dtype=np.int64), {}

        selected_targetids = targetids[rows]
        spectra = {}
        for camera in ("b", "r", "z"):
            if camera_names is not None and camera not in camera_names:
                continue

            prefix = camera.upper()
            wave_ext = f"{prefix}_WAVELENGTH"
            flux_ext = f"{prefix}_FLUX"
            ivar_ext = f"{prefix}_IVAR"
            mask_ext = f"{prefix}_MASK"
            if wave_ext not in extnames or flux_ext not in extnames:
                continue

            spectra[camera] = {
                "wave": hdus[wave_ext].read(),
                "flux": _read_image_rows(hdus[flux_ext], rows),
                "ivar": _read_image_rows(hdus[ivar_ext], rows)
                if ivar_ext in extnames
                else None,
                "mask": _read_image_rows(hdus[mask_ext], rows)
                if mask_ext in extnames
                else None,
            }

    return selected_targetids, spectra


def stack_spectra(
    spectra_paths: list[str],
    redshifts: dict[int, float],
    rest_wave_grid: np.ndarray,
    cameras: list[str] | None,
) -> tuple[np.ndarray, np.ndarray, int]:
    flux_sum = np.zeros(rest_wave_grid.size, dtype=float)
    flux_count = np.zeros(rest_wave_grid.size, dtype=np.int32)
    nspec_used = 0
    wanted_targetids = set(redshifts)

    for spectra_path in spectra_paths:
        targetids, spectra = read_desi_spectra_rows(
            spectra_path, wanted_targetids=wanted_targetids, cameras=cameras
        )
        if len(targetids) == 0 or not spectra:
            continue

        for ispec, targetid in enumerate(targetids):
            if int(targetid) not in redshifts:
                continue

            spectrum_flux_sum = np.zeros(rest_wave_grid.size, dtype=float)
            spectrum_flux_count = np.zeros(rest_wave_grid.size, dtype=np.int16)
            redshift = redshifts[int(targetid)]
            for camera_data in spectra.values():
                ivar = camera_data["ivar"]
                mask = camera_data["mask"]
                add_camera_to_spectrum_grid(
                    rest_wave_grid=rest_wave_grid,
                    obs_wave=camera_data["wave"],
                    flux=camera_data["flux"][ispec],
                    redshift=redshift,
                    spectrum_flux_sum=spectrum_flux_sum,
                    spectrum_flux_count=spectrum_flux_count,
                    ivar=ivar[ispec] if ivar is not None else None,
                    mask=mask[ispec] if mask is not None else None,
                )

            good = spectrum_flux_count > 0
            if np.any(good):
                flux_sum[good] += spectrum_flux_sum[good] / spectrum_flux_count[good]
                flux_count[good] += 1
                nspec_used += 1

    mean_flux = np.full(rest_wave_grid.size, np.nan, dtype=float)
    good = flux_count > 0
    mean_flux[good] = flux_sum[good] / flux_count[good]
    return mean_flux, flux_count, nspec_used


def write_stack(
    output_path: str | Path,
    rest_wave_grid: np.ndarray,
    mean_flux: np.ndarray,
    flux_count: np.ndarray,
    nspec_used: int,
    min_contrib: int,
    clobber: bool,
) -> None:
    output_path = Path(output_path)
    output_flux = mean_flux.copy()
    output_flux[flux_count < min_contrib] = np.nan

    columns = fits.ColDefs(
        [
            fits.Column(
                name="WAVELENGTH", format="D", unit="Angstrom", array=rest_wave_grid
            ),
            fits.Column(name="FLUX", format="D", array=output_flux),
            fits.Column(name="NCONTRIB", format="J", array=flux_count),
        ]
    )
    hdu = fits.BinTableHDU.from_columns(columns, name="STACK")
    hdu.header["NSPEC"] = nspec_used
    hdu.header["MINCONT"] = min_contrib
    hdu.header["STACK"] = "MEAN"
    hdu.header["FRAME"] = "REST"

    hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdul.writeto(output_path, overwrite=clobber)


def main() -> None:
    args = parse_args()
    if args.min_contrib < 1:
        raise ValueError("--min-contrib must be >= 1")

    redshifts = read_redshifts(args.redshifts, args.targetids)
    if not redshifts:
        raise RuntimeError("no usable redshifts found")

    rest_wave_grid = make_wave_grid(
        args.rest_wave_min, args.rest_wave_max, args.rest_wave_step
    )
    mean_flux, flux_count, nspec_used = stack_spectra(
        spectra_paths=args.spectra,
        redshifts=redshifts,
        rest_wave_grid=rest_wave_grid,
        cameras=args.cameras,
    )
    if nspec_used == 0:
        raise RuntimeError("no spectra matched the redshift table and camera selection")

    write_stack(
        output_path=args.output,
        rest_wave_grid=rest_wave_grid,
        mean_flux=mean_flux,
        flux_count=flux_count,
        nspec_used=nspec_used,
        min_contrib=args.min_contrib,
        clobber=args.clobber,
    )


if __name__ == "__main__":
    main()
