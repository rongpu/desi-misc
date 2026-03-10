#!/usr/bin/env python

# modified from https://github.com/araichoor/desihiz/blob/main/bin/desi_custom_coadds

import os
from glob import glob
from datetime import datetime
import numpy as np
import fitsio
from astropy.io import fits
from astropy.table import Table, vstack, join
from pathlib import Path
from desispec.io import findfile
from desimodel.footprint import radec2pix
from desispec.pixgroup import get_exp2healpix_map
from desitarget.targetmask import desi_mask
from desitarget.geomask import match_to
from desisurvey.utils import yesno
from desiutil.log import get_logger
from argparse import ArgumentParser

allowed_steps = ["exptab", "clean_for_copyprod", "copyprod", "prepare"]
nest = True
log = get_logger()


def parse():
    parser = ArgumentParser()
    parser.add_argument("--specprod", help="input specprod (default=daily)", type=str, default="daily")
    parser.add_argument("--my_specprod", help="output specprod (default=daily)", type=str, default="daily")
    parser.add_argument("--my_spectro_redux", help="output folder (default=$DESI_ROOT/users/raichoor/laelbg)", type=str, default=None) 
    parser.add_argument("--my_hpxdir", help="folder where processed files will be: {outdir}/{my_specprod}/healpix/{my_hpxdir} (default=tertiary{prognum}-thru{lastnight})", type=str, default=None)
    parser.add_argument("--my_steps", help="comma-separated list from: {} (default=exptab,copyprod,prepare)".format(",".join(allowed_steps)), type=str, default="exptab,copyprod,prepare")
    parser.add_argument("--prognum", help="tertiary prognum; at least --prognum or --tileids should be set", type=int, default=None)
    parser.add_argument("--tileids", help="comma-separated list of TILEIDs to reduce; at least --prognum or --tileids should be set", type=str, default=None)
    parser.add_argument("--lastnight", help="lastnight", type=int, default=None)
    parser.add_argument("--my_nside", help="healpix nside; has to be >=64 (default=64)", type=int, default=64)
    parser.add_argument("--black_expids", help="comma-separated list of expids to discard", type=str, default=None)
    parser.add_argument("--desispec_version", help="desispec version to load (default=0.67.0, used for loa)", type=str, default="0.67.0")
    parser.add_argument("--overwrite", help="re-generate existing files?", action="store_true")
    parser.add_argument("--cameras", help="comma-separated list of cameras (default=b,r,z)", type=str, default="b,r,z")
    args = parser.parse_args()
    if args.my_spectro_redux is None:
        args.my_spectro_redux = os.path.join(os.getenv("DESI_ROOT"), "users", "raichoor", "laelbg")
    if args.my_hpxdir is None:
        args.my_hpxdir = "tertiary{}-thru{}".format(args.prognum, args.lastnight)
    assert(np.all(np.in1d(args.my_steps.split(","), allowed_steps)))
    assert args.my_nside >= 64
    for kwargs in args._get_kwargs():
        log.info("{}:\t{}".format(kwargs[0], kwargs[1]))
    if (args.prognum is None) and (args.tileids is None):
        msg = "at least one of --prognum or --tileids has to be set"
        log.error(msg)
        raise ValueError(msg)
    return args


def get_exptab_fn(my_spectro_redux, my_specprod, my_hpxdir):
    return os.path.join(my_spectro_redux, my_specprod, "healpix", my_hpxdir, "exposures.fits")


def get_prod_expfn(specprod):
    fn = os.path.join(os.getenv("DESI_ROOT"), "spectro", "redux", specprod, "exposures-{}.fits".format(specprod))
    # if no exposures file, use the daily one                                                                                                                                                      
    if not os.path.isfile(fn):
        print("no {}".format(fn)) 
        fn = fn.replace(specprod, "daily")
        print("use {}".format(fn))
    return fn

def create_exptab(specprod, prognum, tileids, lastnight, my_spectro_redux, my_specprod, my_hpxdir, black_expids=None):
    fn = get_prod_expfn(specprod)
    h = fits.open(fn)
    # exposures
    #
    sel = h["EXPOSURES"].data["NIGHT"] <= lastnight
    if prognum is not None:
        sel &= h["EXPOSURES"].data["FAPRGRM"] == "tertiary{}".format(prognum)
    #
    if tileids is not None:
        sel &= np.in1d(h["EXPOSURES"].data["TILEID"], [int(_) for _ in tileids.split(",")])
    #
    if black_expids is not None:
        sel &= ~np.in1d(h["EXPOSURES"].data["EXPID"], [int(_) for _ in black_expids.split(",")])
    h["EXPOSURES"].data = h["EXPOSURES"].data[sel]
    #
    # frames
    sel = np.in1d(h["FRAMES"].data["EXPID"], h["EXPOSURES"].data["EXPID"])
    h["FRAMES"].data = h["FRAMES"].data[sel]
    # write fits, ascii
    fn = get_exptab_fn(my_spectro_redux, my_specprod, my_hpxdir)
    log.info("writing {}".format(fn))
    h.writeto(fn)
    d = Table.read(fn, "EXPOSURES")
    d.write(fn.replace(".fits", ".ascii"), format="ascii.commented_header")
    log.info("black_expids:\t{}".format(black_expids))
    log.info("selected {} expids observed from {} to {}".format(len(d), d["NIGHT"].min(), d["NIGHT"].max()))
    log.info("selected expids:\t{}".format(d["EXPID"]))


def run_clean_for_copyprod(specprod, my_spectro_redux, my_specprod, my_hpxdir):
    inspecdir = os.path.join(os.getenv("DESI_ROOT"), "spectro", "redux", specprod)
    exptab_fn = get_exptab_fn(my_spectro_redux, my_specprod, my_hpxdir)
    outspecdir = os.path.join(my_spectro_redux, my_specprod)
    # AR list (night, expid)
    d = Table.read(exptab_fn, "EXPOSURES")
    linkfns = []
    # AR calibnight: all files for the night
    # AR exposure_table: ignore, as those are not links
    for night in np.unique(d["NIGHT"]):
        linkfns += sorted(glob(os.path.join(outspecdir, "calibnight", str(night), "*")))
        # linkfns += [os.path.join(outspecdir, "exposure_tables", str(night // 100), "exposure_table_{}.csv".format(night))]
    # AR preproc, exposures: all the files for the expid
    for night, expid in zip(d["NIGHT"], d["EXPID"]):
        linkfns += sorted(glob(os.path.join(outspecdir, "preproc", str(night), "{:08d}".format(expid), "*")))
        linkfns += sorted(glob(os.path.join(outspecdir, "exposures", str(night), "{:08d}".format(expid), "*")))
    #
    log.info("{} symlinks to delete".format(len(linkfns)))
    for linkfn in linkfns:
        fn_mtime = datetime.fromtimestamp(os.lstat(linkfn).st_mtime)
        src_fn = os.readlink(linkfn)
        log.info("{}\t{}".format(fn_mtime, src_fn))
    if yesno("\nok to delete the above {} symlinks?\n".format(len(linkfns))):
        log.info("delete things then!")
        for linkfn in linkfns:
            os.remove(linkfn)
        log.info("deleting done.")
    else:
        log.info("not deleting symlinks, as answer is 'n'")


def run_copyprod(specprod, my_spectro_redux, my_specprod, my_hpxdir):
    inspecdir = os.path.join(os.getenv("DESI_ROOT"), "spectro", "redux", specprod)
    exptab_fn = get_exptab_fn(my_spectro_redux, my_specprod, my_hpxdir).replace(".fits", ".ascii")
    outspecdir = os.path.join(my_spectro_redux, my_specprod)
    cmd = "copyprod --explist {} --abspath {} {}".format(exptab_fn, inspecdir, outspecdir)
    log.info("running: {}".format(cmd))
    _ = os.system(cmd)


def get_children_pixs(pix, nside_in, nside_out):
    npix_ratio = nside_out * nside_out // nside_in // nside_in
    pixlo, pixhi = pix * npix_ratio, (pix+1) * npix_ratio
    pixouts = np.arange(pixlo, pixhi, dtype=int)
    return pixouts


def get_myfns(outdir, outroot, cframe_fns):
    myfns = {
        "cframes"   : " ".join(cframe_fns),
        "spectra"   : os.path.join(outdir, "spectra-{}.fits".format(outroot)),
        "spectralog": os.path.join(outdir, "logs", "spectra-{}.log".format(outroot)),
        "coadd"     : os.path.join(outdir, "coadd-{}.fits".format(outroot)),
        "coaddlog"  : os.path.join(outdir, "logs", "coadd-{}.log".format(outroot)),
        "redrock"   : os.path.join(outdir, "redrock-{}.fits".format(outroot)),
        "rrdetails" : os.path.join(outdir, "rrdetails-{}.h5".format(outroot)),
        "rrlog"     : os.path.join(outdir, "logs", "redrock-{}.log".format(outroot)),
        "qsomgii"   : os.path.join(outdir, "qso_mgii-{}.fits".format(outroot)),
        "qsomgiilog": os.path.join(outdir, "logs", "qso_mgii-{}.log".format(outroot)),
        "qsoqn"     : os.path.join(outdir, "qso_qn-{}.fits".format(outroot)),
        "qsoqnlog"  : os.path.join(outdir, "logs", "qso_qn-{}.log".format(outroot)),
        "emline"    : os.path.join(outdir, "emline-{}.fits".format(outroot)),
        "emlinelog" : os.path.join(outdir, "logs", "emline-{}.log".format(outroot)),
    }
    return myfns


def get_cframe(expid, camera, petal, night):
    fn = findfile("cframe", night=night, expid=expid, camera=camera+str(petal))
    if not os.path.isfile(fn):
        log.warning("missing {}".format(fn))
    return fn


def custom_write(myfs, myfn, step, cmd, overwrite=False):
    myfs[step]["all"].write(cmd)
    #
    if os.path.isfile(os.path.expandvars(myfn)) and (not overwrite):
        log.info("{} already exists".format(myfn))
        pass
    else:
        myfs[step]["rerun"].write(cmd)


# AR loop on healpix pixels
def write_per_step_files(d, cameras, survey, program, myfs, my_spectro_redux, my_specprod, my_hpxdir, my_nside, overwrite=False):
    for i in range(len(d)):
        # 
        healpix = d["HEALPIX"][i]
        expids = [int(expid) for expid in d["EXPIDS"][i].split(",")]
        nights = [int(night) for night in d["NIGHTS"][i].split(",")]
        spectros = [int(spectro) for spectro in d["SPECTROS"][i].split(",")]
        outdir_i = os.path.join(my_spectro_redux, my_specprod, "healpix", my_hpxdir)
        outroot_i = "{}".format(healpix)
        log.info("outdir_i={}, outroot_i={}".format(outdir_i, outroot_i))
        for mydir in [outdir_i, os.path.join(outdir_i, "logs")]:
            if not os.path.isdir(mydir):
                log.info("creating {}".format(mydir))
                Path(os.path.expandvars(mydir)).mkdir(parents=True, exist_ok=True)
        # ====================== cframes ==========================================
        cframe_fns = []
        nexp = 0
        for expid, night, spectro in zip(expids, nights, spectros):
            expid_cframe_fns = np.array([get_cframe(expid, camera, spectro, night) for camera in cameras])
            expid_cframe_isfiles = np.array([os.path.isfile(fn) for fn in expid_cframe_fns])
            if expid_cframe_isfiles.sum() == len(cameras):
                nexp += 1
                cframe_fns += expid_cframe_fns.tolist()
        log.info(
            "{}/{} exposures have brz cframes for healpix={}".format(
                nexp, len(expids), healpix,
            )
        )
        # ======================  myfns        ===================================
        myfns = get_myfns(outdir_i, outroot_i, cframe_fns)
        # ====================== group_spectra ====================================
        cmd = "srun -N 1 -n 1 -c 64 desi_group_spectra --inframes {} --outfile {}".format(
            myfns["cframes"], myfns["spectra"],
        )
        cmd += " --nside {} --healpix {}".format(my_nside, healpix)
        cmd += " --header SPGRP={} SPGRPVAL={} HPXPIXEL={} HPXNSIDE={} HPXNEST={} SURVEY={} PROGRAM={}".format("healpix", healpix, healpix, my_nside, True, survey, program)
        cmd += " &> {}\n".format(myfns["spectralog"])
        custom_write(myfs, myfns["spectra"], "spectra", cmd, overwrite)
        # ====================== coadd_spectra ====================================
        cmd = "desi_coadd_spectra -i {} -o {} &> {}\n".format(
            myfns["spectra"], myfns["coadd"], myfns["coaddlog"],
        )
        custom_write(myfs, myfns["coadd"], "coadd", cmd, overwrite)
        # ====================== redrock ==========================================
        tmax = "01:00:00"
        # cmd = "srun -N 1 -n 64 -c 2 -t {} rrdesi_mpi -i {} -o {} --details {} &> {}\n".format(
        #     tmax, myfns["coadd"], myfns["redrock"], myfns["rrdetails"], myfns["rrlog"],
        cmd = "srun -N 1 -n 64 -c 2 -t {} rrdesi_mpi -i {} -o {} &> {}\n".format(
            tmax, myfns["coadd"], myfns["redrock"], myfns["rrlog"],
        )
        custom_write(myfs, myfns["redrock"], "redrock", cmd, overwrite)
        # ====================== qso_mgii =========================================
        cmd = "srun -N 1 -n 1 -c 128 desi_qso_mgii_afterburner --coadd {} --redrock {} --output {} --target_selection all_targets --save_target all &> {}\n".format(
            myfns["coadd"], myfns["redrock"], myfns["qsomgii"], myfns["qsomgiilog"],
        )
        custom_write(myfs, myfns["qsomgii"], "qsomgii", cmd, overwrite)
        # ====================== qso_qn ===========================================
        cmd = "srun -N 1 -n 1 -c 128 desi_qso_qn_afterburner --coadd {} --redrock {} --output {} --target_selection all_targets --save_target all &> {}\n".format(
            myfns["coadd"], myfns["redrock"], myfns["qsoqn"], myfns["qsoqnlog"],
        )
        custom_write(myfs, myfns["qsoqn"], "qsoqn", cmd, overwrite)
        # ====================== emlinefit ========================================
        cmd = "desi_emlinefit_afterburner --coadd {} --redrock {} --output {} &> {}\n".format(
            myfns["coadd"], myfns["redrock"], myfns["emline"], myfns["emlinelog"],
        )
        custom_write(myfs, myfns["emline"], "emline", cmd, overwrite)


# merging all info in one cat, cutting on secondaries-only
# OUTDIR="$1" TERTIARY="$2" LASTNIGHT="$3" ROOT="$4" python - <<END
def mergecat(mydir, outbasefn, myroot=""):
    # get redrock files + stack
    fns = sorted(glob(os.path.join(mydir, "redrock-*{}.fits".format(myroot))))
    log.info("start working from:\n{}".format("\n".join(fn)))
    d  = vstack([Table(fitsio.read(fn, "REDSHIFTS")) for fn in fns])
    # add fibermap
    fm = vstack([Table(fitsio.read(fn, "FIBERMAP")) for fn in fns])
    d  = join(d, fm, keys=["TARGETID"], metadata_conflicts="silent")
    # add tsnr2
    ts = vstack([Table(fitsio.read(fn, "TSNR2")) for fn in fns])
    d  = join(d, ts, keys=["TARGETID"], metadata_conflicts="silent")
    # add emline
    em = vstack([Table(fitsio.read(fn.replace("redrock", "emline"), "EMLINEFIT")) for fn in fns])
    keys = [key for key in em.dtype.names if key in d.dtype.names and key != "TARGETID"]
    em.remove_columns(keys)
    d  = join(d, em, keys=["TARGETID"], metadata_conflicts="silent")
    # add the GOOD_{BGS,LRG,ELG,QSO} columns
    proceed = True
    for fn in fns:
        for name in ["qso_mgii", "qso_qn", "emline"]:
            name_fn = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace("redrock-", "{}-".format(name)))
            if not os.path.isfile(name_fn):
                proceed = False
    if proceed:
        from desispec.validredshifts import validate as zspec_validate
        valid_d = vstack([zspec_validate(fn, fiberstatus_cut=False) for fn in fns])
        ii = match_to(valid_d["TARGETID"], d["TARGETID"])
        assert(ii.size == len(d))
        assert(np.all(valid_d["TARGETID"][ii] == d["TARGETID"]))
        for key in ["GOOD_BGS", "GOOD_LRG", "GOOD_ELG", "GOOD_QSO"]:
            d[key] = valid_d[key][ii]
    # cut on secondary-only
    sel = d["DESI_TARGET"] == desi_mask["SCND_ANY"]
    log.info("keep {}/{} SCND_ANY targets".format(sel.sum(), len(d)))
    d = d[sel]
    # write
    d.write(os.path.join(mydir, outbasefn))
    

# AR bash wrapper script
def write_bash_wrapper(outfn, my_spectro_redux, my_specprod, myfs, desispec_version, steps, nodes=4):
    f = open(outfn, "w")
    f.write("#!/bin/bash\n")
    f.write("\n") 
    f.write("# Initial setup:\n")
    f.write("source $CFS/desi/software/desi_environment.sh main\n")
    f.write("module swap desispec/{}\n".format(desispec_version))
    # f.write("export PATH=$HOME/software_dev/desispec_healpix/bin:$PATH\n")
    # f.write("export PYTHONPATH=$HOME/software_dev/desispec_healpix/py:$PYTHONPATH\n")
    f.write("export DESI_SPECTRO_REDUX={}\n".format(my_spectro_redux))
    f.write("export SPECPROD={}\n".format(my_specprod))
    f.write("\n")
    f.write("module load parallel\n")
    f.write("salloc --nodes {} --qos interactive --time 4:00:00 --constraint cpu\n".format(nodes))
    f.write("\n")
    for step in steps:
        f.write("parallel --jobs {} --delay 1 < desi_{}.txt; ".format(nodes, step))
    f.write("exit\n")
    f.write("\n")
    f.close()


def main():

    args = parse()

    log.info("output folder: {}".format(os.path.join(args.my_spectro_redux, args.my_specprod, "healpix", args.my_hpxdir)))

    # AR exptab file
    if "exptab" in args.my_steps.split(","):
        create_exptab(args.specprod, args.prognum, args.tileids, args.lastnight, args.my_spectro_redux, args.my_specprod, args.my_hpxdir, black_expids=args.black_expids)

    # AR clean_for_copyprod
    # AR step to remove exisitng symlinks, if we plan to run copyprod after
    # AR typical use:
    # AR    - in my first usage, I was running copyprod without the --abspath argument
    # AR        hence links were relative paths
    # AR    - but in the meantime, I ve moved the folder place
    # AR    - hence the links are broken
    # AR    - so if I want to rerun another reduction, I need to first delete existing links
    if "clean_for_copyprod" in args.my_steps.split(","):
        run_clean_for_copyprod(args.specprod, args.my_spectro_redux, args.my_specprod, args.my_hpxdir)

    # AR copyprod
    if "copyprod" in args.my_steps.split(","):
        run_copyprod(args.specprod, args.my_spectro_redux, args.my_specprod, args.my_hpxdir)

    # AR settings
    if "prepare" in args.my_steps.split(","):
        steps = ["spectra", "coadd", "redrock", "qsomgii", "qsoqn", "emline"]
        os.environ["DESI_SPECTRO_REDUX"] = args.my_spectro_redux
        os.environ["SPECPROD"] = args.my_specprod
        # survey, program = "special", "other"
        exptab_fn = get_exptab_fn(args.my_spectro_redux, args.my_specprod, args.my_hpxdir)
        d = Table.read(exptab_fn, "EXPOSURES")
        # identify the pixels with data
        if args.my_nside != 64:
            tileids = np.unique(d["TILEID"])
            ras, decs = [], []
            # any observed tile has its fiberassign file svn-committed
            for tileid in tileids:
                tileidpad = "{:06d}".format(tileid)
                fn = os.path.join(
                    os.getenv("DESI_TARGET"),
                    "fiberassign",
                    "tiles",
                    "trunk",
                    tileidpad[:3],
                    "fiberassign-{}.fits.gz".format(tileidpad),
                )
                print(fn)
                _ = fitsio.read(fn, ext="FIBERASSIGN", columns=["TARGET_RA", "TARGET_DEC"])
                ras += _["TARGET_RA"].tolist()
                decs += _["TARGET_DEC"].tolist()
            data_unq_pixs = np.unique(radec2pix(args.my_nside, np.array(ras), np.array(decs)))
        #
        assert np.unique(d["SURVEY"]).size == 1
        assert np.unique(d["PROGRAM"]).size == 1
        survey, program = np.unique(d["SURVEY"])[0], np.unique(d["PROGRAM"])[0]
        d = get_exp2healpix_map(exptab_fn, survey=survey, program=program, specprod_dir=os.path.join(args.my_spectro_redux, args.my_specprod))
        myd = {key : [] for key in ["HEALPIX", "EXPIDS", "NIGHTS", "SPECTROS"]}
        nside_in = 64
        for unq_pix in np.unique(d["HEALPIX"]):
            sel = d["HEALPIX"] == unq_pix
            if args.my_nside == 64:
                pixs = [unq_pix]
            else:
                pixs = get_children_pixs(unq_pix, nside_in, args.my_nside)
                pixs = pixs[np.in1d(pixs, data_unq_pixs)]
            for pix in pixs:
                myd["HEALPIX"] += [pix]
                myd["EXPIDS"].append(",".join(d["EXPID"][sel].astype(str)))
                myd["NIGHTS"].append(",".join(d["NIGHT"][sel].astype(str)))
                myd["SPECTROS"].append(",".join(d["SPECTRO"][sel].astype(str)))
        d = Table()
        for key in myd:
            d[key] = myd[key]
        print(d["HEALPIX", "EXPIDS"])

        fn = get_prod_expfn(args.specprod)
        exposures = Table.read(fn, "EXPOSURES")

        # AR open files for writing
        myfs = {}
        for step in steps:
            myfs[step] = {
                "all" : open("desi_{}_all.txt".format(step), "w"),
                "rerun" : open("desi_{}.txt".format(step), "w"),
            }

        # AR write files
        write_per_step_files(d, args.cameras.split(","), survey, program, myfs, args.my_spectro_redux, args.my_specprod, args.my_hpxdir, args.my_nside, overwrite=args.overwrite)

        # AR close files
        for step in steps:
            myfs[step]["all"].close()
            myfs[step]["rerun"].close()

        # AR bash wrapper script
        bashfn = "laelbg-custom_coadds_launch.sh"
        write_bash_wrapper(bashfn, args.my_spectro_redux, args.my_specprod, myfs, args.desispec_version, steps, nodes=4)
        os.system("cat {}".format(bashfn))

if __name__ == "__main__":
    main()

