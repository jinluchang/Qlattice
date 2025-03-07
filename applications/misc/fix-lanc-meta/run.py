#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import cevec_io_meta

import os

from jobs import *

load_path_list[:] = [
        "results",
        "../mk-gf-gt/results",
        os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-lanc-no-meta/results"),
        "/sdcc/u/jluchang/qcdqedta/luchang/data-gen/lanc/32Dfine-metadata/results-no-meta/",
        ]

def prod(l):
    p = 1
    for x in l:
        p *= x
    return p

@q.timer_verbose
def load_eig(path, job_tag, inv_type = 0, inv_acc = 0):
    total_site = ru.get_total_site(job_tag)
    fermion_params = ru.get_param_fermion(job_tag, inv_type, inv_acc)
    grids = qg.get_fgrid(total_site, fermion_params)
    eig = cevec_io_meta.load(path, grids = grids)
    return eig

@q.timer_verbose
def save_ceig(path, eig, job_tag, inv_type = 0, inv_acc = 0, *, crc32 = None, mpi = None):
    if path is None:
        return
    save_params = ru.get_clanc_params(job_tag, inv_type, inv_acc)["save_params"]
    nsingle = save_params["nsingle"]
    if mpi is None:
        mpi = save_params["mpi"]
    fmt = g.format.cevec({"nsingle": nsingle, "mpi": [ 1 ] + mpi, "max_read_blocks": 8})
    cevec_io_meta.save_meta(path, eig, fmt.params, crc32 = crc32);

def save_metadata(path, params, inv_type, inv_acc, *, mpi = None):
    filename = path
    fermion_params = params["fermion_params"][inv_type][inv_acc]
    ls = ru.get_ls_from_fermion_params(fermion_params)
    total_site = params["total_site"]
    cparams = params["clanc_params"][inv_type][inv_acc]
    neigen = cparams["irl_params"]["Nstop"]
    nbasis = cparams["nbasis"]
    nsingle = cparams["save_params"]["nsingle"]
    if mpi is None:
        mpi = cparams["save_params"]["mpi"]
    block = cparams["block"]
    s = [ total_site[i] // mpi[i] for i in range(4) ] + [ ls, ]
    b = block + [ ls, ]
    nb = [ s[i] // b[i] for i in range(5) ]
    nd = 5
    nrank = prod(mpi) 
    blocks = prod(nb)
    crc32 = [ 0, ] * nrank
    FP16_COEF_EXP_SHARE_FLOATS = 10
    # write metadata
    if q.get_id_node() == 0:
        fmeta = open("%s/metadata.txt" % filename, "wt")
        for i in range(nd):
            fmeta.write("s[%d] = %d\n" % (i, s[i]))
        for i in range(nd):
            fmeta.write("b[%d] = %d\n" % (i, b[i]))
        for i in range(nd):
            fmeta.write("nb[%d] = %d\n" % (i, nb[i]))
        fmeta.write("neig = %d\n" % neigen)
        fmeta.write("nkeep = %d\n" % nbasis)
        fmeta.write("nkeep_single = %d\n" % nsingle)
        fmeta.write("blocks = %d\n" % blocks)
        fmeta.write("FP16_COEF_EXP_SHARE_FLOATS = %d\n" % FP16_COEF_EXP_SHARE_FLOATS)
        for i in range(len(crc32)):
            fmeta.write("crc32[%d] = %X\n" % (i, crc32[i]))
        fmeta.close()

@q.timer_verbose
def run_eig_fix_meta(job_tag, traj, get_gf, inv_type = 0, inv_acc = 0, *, mpi_original = None):
    assert get_gf is not None
    gf = get_gf()
    path_eig = get_load_path(f"{job_tag}/eig/traj-{traj}")
    assert path_eig is not None
    save_metadata(path_eig, rup.dict_params[job_tag], inv_type, inv_acc, mpi = mpi_original)
    basis, cevec, crc32 = load_eig(path_eig, job_tag, inv_type, inv_acc)
    for i in range(len(crc32)):
        crc32[i] = q.glb_sum(crc32[i])
    smoothed_evals = ru.get_smoothed_evals(basis, cevec, gf, job_tag, inv_type, inv_acc)
    q.displayln_info("smoothed_evals=", smoothed_evals)
    eig = basis, cevec, smoothed_evals
    save_ceig(path_eig, eig, job_tag, inv_type, inv_acc, crc32 = crc32, mpi = mpi_original)
    test_eig(gf, eig, job_tag, inv_type)

@q.timer_verbose
def run_eig_fix_reshape(job_tag, traj, get_gf, inv_type = 0, inv_acc = 0, *, mpi_original = None):
    assert get_gf is not None
    gf = get_gf()
    path = f"{job_tag}/eig/traj-{traj}"
    path_eig = get_load_path(path)
    assert path_eig is not None
    save_metadata(path_eig, rup.dict_params[job_tag], inv_type, inv_acc, mpi = mpi_original)
    basis, cevec, crc32 = load_eig(path_eig, job_tag, inv_type, inv_acc)
    smoothed_evals = ru.get_smoothed_evals(basis, cevec, gf, job_tag, inv_type, inv_acc)
    q.displayln_info("smoothed_evals=", smoothed_evals)
    eig = basis, cevec, smoothed_evals
    ru.save_ceig(get_save_path(path + ".partial"), eig, job_tag, inv_type, inv_acc);
    q.qrename_info(get_save_path(path + ".partial"), get_save_path(path))
    test_eig(gf, eig, job_tag, inv_type)

def guess_eig_mpi(job_tag, traj):
    if job_tag != "32Dfine":
        mpi = rup.dict_params[job_tag]["lanc-mpi-original"]
    else:
        assert job_tag == "32Dfine"
        path_eig = get_load_path(f"{job_tag}/eig/traj-{traj}")
        if path_eig is None:
            mpi = rup.dict_params[job_tag]["lanc-mpi-original"]
        elif q.does_file_exist_sync_node(os.path.join(path_eig, "31/0000000511.compressed")):
            num_node = 512
            mpi = [ 4, 4, 4, 8, ]
        elif q.does_file_exist_sync_node(os.path.join(path_eig, "31/0000000255.compressed")):
            num_node = 256
            mpi = [ 4, 4, 4, 4, ]
        elif q.does_file_exist_sync_node(os.path.join(path_eig, "31/0000000127.compressed")):
            num_node = 128
            mpi = [ 4, 4, 4, 2, ]
        elif q.does_file_exist_sync_node(os.path.join(path_eig, "31/0000000063.compressed")):
            num_node = 64
            mpi = [ 2, 2, 2, 8, ]
        else:
            assert False
    q.displayln_info(f"guess_eig_mpi: {job_tag} {traj} mpi={mpi}")
    return mpi

def run_eig_fix(job_tag, traj, get_gf, inv_type = 0, inv_acc = 0):
    if q.obtain_lock(f"locks/{job_tag}-{traj}-run-eig-fix"):
        run_eig_fix_meta(job_tag, traj, get_gf, inv_type, inv_acc, mpi_original = guess_eig_mpi(job_tag, traj))
        # run_eig_fix_reshape(job_tag, traj, get_gf, inv_type, inv_acc, mpi_original = guess_eig_mpi(job_tag, traj))
        q.release_lock()

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/eig/traj-{traj}/eigen-values.txt",
            f"{job_tag}/eig/traj-{traj}/metadata.txt",
            ]
    fns_need = [
            (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/eig/traj-{traj}",
            ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_gf = run_gf(job_tag, traj)
    #
    run_eig_fix(job_tag, traj, get_gf)
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["32Dfine"]["trajs"] = list(range(200, 5000, 10))

rup.dict_params["test-4nt8"]["lanc-mpi-original"] = [ 1, 1, 1, 4, ]
rup.dict_params["test-4nt16"]["lanc-mpi-original"] = [ 1, 1, 1, 4, ]
rup.dict_params["32Dfine"]["lanc-mpi-original"] = [ 2, 4, 4, 4, ]

# rup.dict_params["32Dfine"]["load_config_params"]["twist_boundary_at_boundary"] = [ 0.0, 0.0, 0.0, 0.0, ]

# rup.dict_params["32Dfine"]["clanc_params"][0][0]["smoother_params"]["maxiter"] = 10

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "test-8nt16",
        # "test-16nt32",
        # "test-32nt64",
        # "test-48nt96",
        # "test-64nt128",
        # "test-96nt192",
        # "test-128nt256",
        # "24D",
        # "32Dfine",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
