#!/usr/bin/env python3

from auto_contractor.distillation_mat_op import *

q.begin_with_mpi()
q.qremove_all_info("results")
q.qmkdir_info("results")

nn_s = 4

nn_c = 16

def make_rnd_sm(rs : q.RngState):
    x = np.zeros((nn_s, nn_s,), dtype = np.complex128)
    rs.u_rand_fill(x, 1.0, -1.0)
    return x

def make_rnd_cm(rs : q.RngState):
    x = np.zeros((nn_c, nn_c,), dtype = np.complex128)
    rs.u_rand_fill(x, 1.0, -1.0)
    return x

def make_rnd_wm(rs : q.RngState):
    x = np.zeros((nn_s, nn_s, nn_c, nn_c,), dtype = np.complex128)
    rs.u_rand_fill(x, 1.0, -1.0)
    return x

sm1 = make_rnd_sm(q.RngState("seed-sm1"))
sm2 = make_rnd_sm(q.RngState("seed-sm2"))
v1 = mat_tr_sm(sm1)
v2 = mat_tr_sm(sm2)
q.displayln_info(f"CHECK: mat_tr_sm: {v1:.10f} {v2:.10f}")

cm1 = make_rnd_cm(q.RngState("seed-cm1"))
cm2 = make_rnd_cm(q.RngState("seed-cm2"))
v1 = mat_tr_cm(cm1)
v2 = mat_tr_cm(cm2)
q.displayln_info(f"CHECK: mat_tr_cm: {v1:.10f} {v2:.10f}")

wm1 = make_rnd_wm(q.RngState("seed-wm1"))
wm2 = make_rnd_wm(q.RngState("seed-wm2"))
v1 = mat_tr_wm(wm1)
v2 = mat_tr_wm(wm2)
q.displayln_info(f"CHECK: mat_tr_wm: {v1:.10f} {v2:.10f}")

v1 = mat_tr_wm_wm(wm1, wm1)
v2 = mat_tr_wm_wm(wm2, wm2)
v3 = mat_tr_wm_wm(wm1, wm2)
q.displayln_info(f"CHECK: mat_tr_wm_wm: {v1:.10f} {v2:.10f} {v3:.10f}")

q.check_all_files_crc32_info("results")
q.timer_display()
q.displayln_info(f"CHECK: finished successfully.")
q.end_with_mpi()
