#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

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
def compute_eig(gf, job_tag, inv_type = 0, inv_acc = 0, *, path = None):
    # return a function ``get_eig''
    # ``get_eig()'' return the ``eig''
    load_eig = ru.load_eig_lazy(get_load_path(path), job_tag)
    if load_eig is not None:
        return load_eig
    # evec, evals = ru.mk_eig(gf, job_tag, inv_type, inv_acc)
    basis, cevec, smoothed_evals = ru.mk_ceig(gf, job_tag, inv_type, inv_acc)
    eig = [ basis, cevec, smoothed_evals ]
    ru.save_ceig(get_save_path(path), eig, job_tag, inv_type, inv_acc);
    def get_eig():
        return eig
    return get_eig

@q.timer
def test_eig(gf, eig, job_tag, inv_type):
    geo = gf.geo()
    src = q.FermionField4d(geo)
    q.displayln_info(f"CHECK: src norm {src.qnorm()}")
    src.set_rand(q.RngState("test_eig:{id(inv)}"))
    sol_ref = ru.get_inv(gf, job_tag, inv_type, inv_acc = 2, eig = eig, eps = 1e-10, mpi_split = False, qtimer = False) * src
    q.displayln_info(f"CHECK: sol_ref norm {sol_ref.qnorm():.14E} with eig")
    for inv_acc in [0, 1, 2]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig, mpi_split = False, qtimer = False) * src
        sol -= sol_ref
        q.displayln_info(f"CHECK: sol diff norm {sol.qnorm():.2E} inv_acc={inv_acc} with eig")
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, mpi_split = False, qtimer = False) * src
        sol -= sol_ref
        q.displayln_info(f"CHECK: sol diff norm {sol.qnorm():.2E} inv_acc={inv_acc} without eig")

@q.timer
def run(job_tag, traj):
    q.check_stop()
    q.check_time_limit()
    #
    q.qmkdir_info(get_save_path(f""))
    q.qmkdir_info(get_save_path(f"eig"))
    q.qmkdir_info(get_save_path(f"eig/{job_tag}"))
    #
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    q.displayln_info("CHECK: geo.show() =", geo.show())
    #
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    if path_gf is None:
        if job_tag[:5] == "test-":
            gf = ru.mk_sample_gauge_field(job_tag, f"{traj}")
            q.qmkdir_info(get_save_path(f"configs"))
            q.qmkdir_info(get_save_path(f"configs/{job_tag}"))
            path_gf = get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
            # gf.save(path_gf)
            qg.save_gauge_field(gf, path_gf)
        else:
            assert False
    gf = ru.load_config(job_tag, path_gf)
    gf.show_info()
    #
    get_eig = compute_eig(gf, job_tag, inv_type = 0, path = f"eig/{job_tag}/traj={traj}")
    test_eig(gf, get_eig(), job_tag, inv_type = 0)

qg.begin_with_gpt()

q.qremove_all_info("results")

job_tag = "test-4nt16"
traj = 1000
q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
q.displayln_info("CHECK: ", rup.dict_params[job_tag])
run(job_tag, traj)

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

qg.end_with_gpt()
