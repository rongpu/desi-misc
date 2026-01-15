#!/usr/bin/env python

import fitsio
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u


# https://github.com/araichoor/desihiz/blob/8ef9a3326de105013d9ca59e57ba988da1b8055b/py/desihiz/hizmerge_io.py#L320
# https://desisurvey.slack.com/archives/C027M1AF31C/p1734727460173929
def get_ext_coeffs():

    return {
        "SYNTHG": 3.240,
        "G": 3.240,
        "R": 2.276,
        "I": 1.633,
        "Z": 1.263,
        "Y": 1.076,
        "M411": 4.290,
        "M438": 4.099,
        "M464": 3.877,
        "M490": 3.634,
        "M517": 3.389,
    }


# [decam-chatter 18711] Suggested MASKBITS for IBIS
# 2024-12-20
def get_maskbits_cut_bits(bands):

    cut_names = ["NPRIMARY", "BRIGHT", "MEDIUM", "GALAXY"]
    for band in bands:
        cut_names += ["SATUR_{}".format(band), "ALLMASK_{}".format(band)]

    # pick one tractor file, to grab header
    fn = "/pscratch/sd/d/dstn/ibis3/tractor/032/tractor-0329m042.fits"
    hdr = fitsio.read_header(fn, 0)
    keys = np.array([key for key in hdr if key[:5] == "MBIT_"])
    names = np.array([hdr[key] for key in keys])

    cut_keys = keys[np.in1d(names, cut_names)]
    cut_bits = np.array([int(key.split("_")[-1]) for key in cut_keys])

    return cut_bits


def main():

    forced = "wide"

    fibmagmin, fibmagmax = 22.0, 25.25
    mb_colmin = 0.3
    ri_colmax = 0.35
    # gr_colmax = 0.7

    outfn = '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/xmm_highz/xmm_lae_placeholder.fits'

    fn = "/dvs_ro/cfs/cdirs/desi/users/raichoor/ibis/ibis3-xmm-forced-hsc-{}.fits".format(forced)

    # AR offset field by +0.3 in dec to avoid shallow m464 region
    field_ra, field_dec, field_rad_deg = 35.75, -4.75, 1.8 - 0.1
    area = np.pi * field_rad_deg ** 2

    fnodet = 0.001
    magnodet = 22.5 - 2.5 * np.log10(fnodet) # 30

    d = Table(fitsio.read(fn))
    # remove ambiguous columns (as those MAG do not have ebv-corr,
    #   and we use different values here; also remove all x-match with desi infos)
    black_keys = []
    black_keys += [_ for _ in d.colnames if _[:4] == "MAG_"]
    black_keys += [_ for _ in d.colnames if _[:7] == "FIBMAG_"]
    black_keys += [_ for _ in d.colnames if "DESI" in _]
    print("remove these columns: {}".format(",".join(black_keys)))
    d.remove_columns(black_keys)

    g_mags = 22.5 - 2.5 * np.log10(np.clip(d["FORCED_FLUX_G"], fnodet, None)) - get_ext_coeffs()["G"] * d["EBV"]
    r_mags = 22.5 - 2.5 * np.log10(np.clip(d["FORCED_FLUX_R"], fnodet, None)) - get_ext_coeffs()["R"] * d["EBV"]
    i_mags = 22.5 - 2.5 * np.log10(np.clip(d["FORCED_FLUX_I"], fnodet, None)) - get_ext_coeffs()["I"] * d["EBV"]

    sel_all = np.full(len(d), False)

    # selband = 'M490'
    # blueband = 'M464'
    for selband, blueband in [['M438', 'M411'], ['M464', 'M438'], ['M490', 'M464'], ['M517', 'M490']]:

        print('\n', selband)

        # cut on radius
        field_c = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame="icrs")                                                                                                   
        cs = SkyCoord(ra=d["RA"] * u.deg, dec=d["DEC"] * u.deg, frame="icrs")
        d["FIELD_RADIUS_OFFSET"] = cs.separation(field_c).deg
        sel = d["FIELD_RADIUS_OFFSET"] < field_rad_deg
        print("{:.1f}k/deg2\tafter radius cut".format(sel.sum() / 1000. / area))

        # cut on maskbits
        cut_bits = get_maskbits_cut_bits(["M464", "M490"])
        for bit in cut_bits:
            sel &= (d["MASKBITS"] & 2**bit) == 0
        print("{:.1f}k/deg2\tafter maskbits cut".format(sel.sum() / 1000. / area))

        # cut on fibmag
        fibmags = {
            band: 22.5 - 2.5 * np.log10(np.clip(d["FIBERFLUX_{}".format(band)], fnodet, None)) - get_ext_coeffs()[band] * d["EBV"]
            for band in [blueband, selband]
        }
        sel &= (fibmags[selband] > fibmagmin) & (fibmags[selband] < fibmagmax)
        print("{:.1f}k/deg2\tafter {} < {} fibmag < {} cut".format(sel.sum() / 1000. / area, fibmagmin, selband, fibmagmax))

        # cut on m464-m490
        sel &= fibmags[blueband] - fibmags[selband] > mb_colmin
        print("{:.1f}k/deg2\tafter {}-{} > {} cut".format(sel.sum() / 1000. / area, blueband, selband, mb_colmin))

        # cut on r-i
        # case: no_r, no_i: will be selected by ri<0.5
        # case: r, no_i: we add it manually here
        # case: no_r, i: t will be selected by ri<0.5; we don t want it, but ok
        sel &= (r_mags - i_mags < ri_colmax) | ((r_mags < magnodet) & (i_mags == magnodet))
        print("{:.1f}k/deg2\tafter r-i < {} cut".format(sel.sum() / 1000. / area, ri_colmax))
        # sel &= (g_mags - r_mags < gr_colmax) | ((g_mags < magnodet) & (r_mags == magnodet))
        # print("{:.1f}k/deg2\tafter g-r < {} cut".format(sel.sum() / 1000. / area, ri_colmax))

        # better safe than sorry... remove a handful of dubious bright r-mag targets
        sel2 = (sel) & (r_mags > 17)
        print("remove {} targets with r_mag < 17".format(((sel) & (~sel2)).sum()))
        sel &= sel2

        d[selband+'_sel'] = sel.copy()
        sel_all |= sel

        print("{:.1f}k/deg2 combined".format(sel_all.sum() / 1000. / area))

    d = d[sel_all]

    # header infos
    d.meta["FIELDRA"], d.meta["FIELDDEC"], d.meta["FIELDRAD"] = field_ra, field_dec, field_rad_deg
    # d.meta["SELBAND"], d.meta["BLUBAND"] = selband, blueband
    # d.meta["BB_COLMIN"], d.meta["RI_COLMAX"] = mb_colmin, ri_colmax

    d.write(outfn, overwrite=True)

if __name__ == "__main__":
    main()                                                                                                                                                                                       

