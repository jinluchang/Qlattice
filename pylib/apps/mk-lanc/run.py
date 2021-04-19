#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru

import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [ "results" ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer
def mk_sample_gauge_field(job_tag, traj):
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_sample_gauge_field")
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    gf = q.GaugeField(geo)
    gf.set_rand(rs, sigma = 0.25, n_step = 4)
    for i in range(4):
        q.gf_wilson_flow_step(gf, 0.05)
    gf.unitarize()
    return gf

@q.timer
def load_eig(job_tag, inv_type, path):
    path_load = get_load_path(path)
    if path_load is None:
        return None
    total_site = ru.get_total_site(job_tag)
    fermion_params = ru.get_fermion_param(job_tag, inv_type, inv_acc = 0)
    eig = g.load(path_load, grids = qg.get_fgrid(total_site, fermion_params))
    return eig

@q.timer
def compute_eig(gf, job_tag, inv_type, *, path = None, nsingle = 10, mpi = [ 1, 1, 1, 4 ]):
    eig = load_eig(job_tag, inv_type, path)
    if eig is not None:
        return eig 
    # evec, evals = ru.mk_eig(gf, job_tag, inv_type)
    basis, cevec, smoothed_evals = ru.mk_ceig(gf, job_tag, inv_type)
    eig = [basis, cevec, smoothed_evals]
    fmt = g.format.cevec({"nsingle": nsingle, "mpi": [ 1 ] + mpi, "max_read_blocks": 8})
    if path is not None:
        g.save(get_save_path(path), eig, fmt);
    return eig

@q.timer
def test_eig(gf, eig, job_tag, inv_type):
    geo = gf.geo()
    src = q.FermionField4d(geo)
    q.displayln_info(f"src norm {src.qnorm()}")
    src.set_rand(q.RngState("test_eig:{id(inv)}"))
    sol_ref = ru.get_inv(gf, job_tag, inv_type, 2, eig = eig, eps = 1e-10) * src
    q.displayln_info(f"sol_ref norm {sol_ref.qnorm()} with eig")
    for inv_acc in [0, 1, 2]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} with eig")
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} without eig")

def run(job_tag, traj):
    q.qmkdir_info(f"locks")
    q.qmkdir_info(get_save_path(f""))
    q.qmkdir_info(get_save_path(f"configs"))
    q.qmkdir_info(get_save_path(f"configs/{job_tag}"))
    q.qmkdir_info(get_save_path(f"eig"))
    q.qmkdir_info(get_save_path(f"eig/{job_tag}"))
    #
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    q.displayln_info("geo.show() =", geo.show())
    #
    gf = mk_sample_gauge_field(job_tag, traj)
    gf.show_info()
    gf.save(get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}"))
    #
    if q.obtain_lock(f"locks/{job_tag}-{traj}-compute-eig"):
        eig = compute_eig(gf, job_tag, inv_type = 0, path = f"eig/{job_tag}/traj={traj}")
        test_eig(gf, eig, job_tag, inv_type = 0)
        q.release_lock()

qg.begin_with_gpt()

for job_tag in [ "test-4nt16" ]:
    for traj in range(1000, 1100, 100):
        run(job_tag, traj)
        q.timer_display()

qg.end_with_gpt()
