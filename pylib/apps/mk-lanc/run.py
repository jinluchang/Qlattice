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
def load_eig_lazy(job_tag, inv_type, path):
    # return ``None'' or a function ``load_eig''
    # ``load_eig()'' return the ``eig''
    path_load = get_load_path(path)
    if path_load is None:
        return None
    total_site = ru.get_total_site(job_tag)
    fermion_params = ru.get_fermion_param(job_tag, inv_type, inv_acc = 0)
    grids = qg.get_fgrid(total_site, fermion_params)
    eig = [ None ]
    def load_eig():
        assert isinstance(eig, list)
        assert len(eig) == 1
        if eig[0] is None:
            eig[0] = g.load(path_load, grids = grids)
        return eig[0]
    return load_eig

@q.timer
def compute_eig(gf, job_tag, inv_type, *, path = None, nsingle = 10, mpi = [ 1, 1, 1, 4 ]):
    # return a function ``get_eig''
    # ``get_eig()'' return the ``eig''
    load_eig = load_eig_lazy(job_tag, inv_type, path)
    if load_eig is not None:
        return load_eig
    # evec, evals = ru.mk_eig(gf, job_tag, inv_type)
    basis, cevec, smoothed_evals = ru.mk_ceig(gf, job_tag, inv_type)
    eig = [ basis, cevec, smoothed_evals ]
    fmt = g.format.cevec({"nsingle": nsingle, "mpi": [ 1 ] + mpi, "max_read_blocks": 8})
    if path is not None:
        g.save(get_save_path(path), eig, fmt);
    def get_eig():
        return eig
    return get_eig

@q.timer
def test_eig(gf, eig, job_tag, inv_type):
    geo = gf.geo()
    src = q.FermionField4d(geo)
    q.displayln_info(f"src norm {src.qnorm()}")
    src.set_rand(q.RngState("test_eig:{id(inv)}"))
    sol_ref = ru.get_inv(gf, job_tag, inv_type, inv_acc = 2, eig = eig, eps = 1e-10, timer = False) * src
    q.displayln_info(f"sol_ref norm {sol_ref.qnorm()} with eig")
    for inv_acc in [0, 1, 2]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} with eig")
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} without eig")

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

def run(job_tag, traj):
    q.qmkdir_info(f"locks")
    q.qmkdir_info(get_save_path(f""))
    q.qmkdir_info(get_save_path(f"eig"))
    q.qmkdir_info(get_save_path(f"eig/{job_tag}"))
    q.qmkdir_info(get_save_path(f"configs"))
    q.qmkdir_info(get_save_path(f"configs/{job_tag}"))
    #
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    q.displayln_info("geo.show() =", geo.show())
    #
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    if path_gf is None:
        gf = mk_sample_gauge_field(job_tag, traj)
        gf.show_info()
        gf.save(get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}"))
    else:
        gf = q.GaugeField()
        gf.load(path_gf)
    #
    if q.obtain_lock(f"locks/{job_tag}-{traj}-compute-eig"):
        get_eig = compute_eig(gf, job_tag, inv_type = 0, path = f"eig/{job_tag}/traj={traj}")
        test_eig(gf, get_eig(), job_tag, inv_type = 0)
        q.release_lock()

qg.begin_with_gpt()

for job_tag in [ "test-4nt16" ]:
    for traj in range(1000, 1400, 100):
        run(job_tag, traj)
        q.timer_display()

qg.end_with_gpt()
