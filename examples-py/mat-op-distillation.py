#!/usr/bin/env python3

import qlat as q
from auto_contractor.runtime_distillation import *

q.begin_with_mpi()
q.qremove_all_info("results")
q.qmkdir_info("results")

nn_s = 4

nn_c = 3

def make_rnd_sm(rs: q.RngState):
    x = np.zeros((nn_s, nn_s,), dtype = np.complex128)
    rs.u_rand_fill(x, 1.0, -1.0)
    return x

def make_rnd_cm(rs: q.RngState):
    x = np.zeros((nn_c, nn_c,), dtype = np.complex128)
    rs.u_rand_fill(x, 1.0, -1.0)
    return x

def make_rnd_wm(rs: q.RngState):
    x = np.zeros((nn_s, nn_c, nn_s, nn_c,), dtype = np.complex128)
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
v4 = mat_tr_wm_wm(wm2, wm1)
q.displayln_info(f"CHECK: mat_tr_wm_wm: {v1:.10f} {v2:.10f} {v3:.10f} {v4:.10f}")

v1 = mat_tr_sm_sm(sm1, sm1)
v2 = mat_tr_sm_sm(sm2, sm2)
v3 = mat_tr_sm_sm(sm1, sm2)
v4 = mat_tr_sm_sm(sm2, sm1)
q.displayln_info(f"CHECK: mat_tr_sm_sm: {v1:.10f} {v2:.10f} {v3:.10f} {v4:.10f}")

v1 = mat_tr_cm_cm(cm1, cm1)
v2 = mat_tr_cm_cm(cm2, cm2)
v3 = mat_tr_cm_cm(cm1, cm2)
v4 = mat_tr_cm_cm(cm2, cm1)
q.displayln_info(f"CHECK: mat_tr_cm_cm: {v1:.10f} {v2:.10f} {v3:.10f} {v4:.10f}")

v1 = mat_tr_wm_sm(wm1, sm1)
q.displayln_info(f"CHECK: mat_tr_wm_sm: {v1:.10f}")

v1 = mat_tr_sm_wm(sm1, wm1)
q.displayln_info(f"CHECK: mat_tr_sm_wm: {v1:.10f}")

v1 = mat_tr_wm_cm(wm1, cm1)
q.displayln_info(f"CHECK: mat_tr_wm_cm: {v1:.10f}")

v1 = mat_tr_cm_wm(cm1, wm1)
q.displayln_info(f"CHECK: mat_tr_cm_wm: {v1:.10f}")

v1 = q.get_double_sig(mat_mul_wm_wm(wm1, wm1), q.RngState("seed-sig-mul"))
v2 = q.get_double_sig(mat_mul_wm_wm(wm2, wm2), q.RngState("seed-sig-mul"))
v3 = q.get_double_sig(mat_mul_wm_wm(wm1, wm2), q.RngState("seed-sig-mul"))
v4 = q.get_double_sig(mat_mul_wm_wm(wm2, wm1), q.RngState("seed-sig-mul"))
q.displayln_info(f"CHECK: mat_mul_wm_wm: {v1:.10f} {v2:.10f} {v3:.10f} {v4:.10f}")
assert abs(mat_tr_wm(mat_mul_wm_wm(wm1, wm1)) - mat_tr_wm_wm(wm1, wm1)) <= 1e-10
assert abs(mat_tr_wm(mat_mul_wm_wm(wm2, wm2)) - mat_tr_wm_wm(wm2, wm2)) <= 1e-10
assert abs(mat_tr_wm(mat_mul_wm_wm(wm1, wm2)) - mat_tr_wm_wm(wm1, wm2)) <= 1e-10
assert abs(mat_tr_wm(mat_mul_wm_wm(wm2, wm1)) - mat_tr_wm_wm(wm2, wm1)) <= 1e-10

v1 = q.get_double_sig(mat_mul_sm_sm(sm1, sm1), q.RngState("seed-sig-mul"))
v2 = q.get_double_sig(mat_mul_sm_sm(sm2, sm2), q.RngState("seed-sig-mul"))
v3 = q.get_double_sig(mat_mul_sm_sm(sm1, sm2), q.RngState("seed-sig-mul"))
v4 = q.get_double_sig(mat_mul_sm_sm(sm2, sm1), q.RngState("seed-sig-mul"))
q.displayln_info(f"CHECK: mat_mul_sm_sm: {v1:.10f} {v2:.10f} {v3:.10f} {v4:.10f}")
assert abs(mat_tr_sm(mat_mul_sm_sm(sm1, sm1)) - mat_tr_sm_sm(sm1, sm1)) <= 1e-10
assert abs(mat_tr_sm(mat_mul_sm_sm(sm2, sm2)) - mat_tr_sm_sm(sm2, sm2)) <= 1e-10
assert abs(mat_tr_sm(mat_mul_sm_sm(sm1, sm2)) - mat_tr_sm_sm(sm1, sm2)) <= 1e-10
assert abs(mat_tr_sm(mat_mul_sm_sm(sm2, sm1)) - mat_tr_sm_sm(sm2, sm1)) <= 1e-10

v1 = q.get_double_sig(mat_mul_cm_cm(cm1, cm1), q.RngState("seed-sig-mul"))
v2 = q.get_double_sig(mat_mul_cm_cm(cm2, cm2), q.RngState("seed-sig-mul"))
v3 = q.get_double_sig(mat_mul_cm_cm(cm1, cm2), q.RngState("seed-sig-mul"))
v4 = q.get_double_sig(mat_mul_cm_cm(cm2, cm1), q.RngState("seed-sig-mul"))
q.displayln_info(f"CHECK: mat_mul_cm_cm: {v1:.10f} {v2:.10f} {v3:.10f} {v4:.10f}")
assert abs(mat_tr_cm(mat_mul_cm_cm(cm1, cm1)) - mat_tr_cm_cm(cm1, cm1)) <= 1e-10
assert abs(mat_tr_cm(mat_mul_cm_cm(cm2, cm2)) - mat_tr_cm_cm(cm2, cm2)) <= 1e-10
assert abs(mat_tr_cm(mat_mul_cm_cm(cm1, cm2)) - mat_tr_cm_cm(cm1, cm2)) <= 1e-10
assert abs(mat_tr_cm(mat_mul_cm_cm(cm2, cm1)) - mat_tr_cm_cm(cm2, cm1)) <= 1e-10

v1 = q.get_double_sig(mat_mul_wm_sm(wm1, sm1), q.RngState("seed-sig-mul"))
q.displayln_info(f"CHECK: mat_mul_wm_sm: {v1:.10f}")
assert abs(mat_tr_wm(mat_mul_wm_sm(wm1, sm1)) - mat_tr_wm_sm(wm1, sm1)) <= 1e-10

v1 = q.get_double_sig(mat_mul_sm_wm(sm1, wm1), q.RngState("seed-sig-mul"))
q.displayln_info(f"CHECK: mat_mul_sm_wm: {v1:.10f}")
assert abs(mat_tr_wm(mat_mul_sm_wm(sm1, wm1)) - mat_tr_sm_wm(sm1, wm1)) <= 1e-10

v1 = q.get_double_sig(mat_mul_wm_cm(wm1, cm1), q.RngState("seed-sig-mul"))
q.displayln_info(f"CHECK: mat_mul_wm_cm: {v1:.10f}")
assert abs(mat_tr_wm(mat_mul_wm_cm(wm1, cm1)) - mat_tr_wm_cm(wm1, cm1)) <= 1e-10

v1 = q.get_double_sig(mat_mul_cm_wm(cm1, wm1), q.RngState("seed-sig-mul"))
q.displayln_info(f"CHECK: mat_mul_cm_wm: {v1:.10f}")
assert abs(mat_tr_wm(mat_mul_cm_wm(cm1, wm1)) - mat_tr_cm_wm(cm1, wm1)) <= 1e-10

v1 = q.get_double_sig(wilson_matrix_g5_herm(wm1), q.RngState("seed-sig-mul"))
v2 = q.get_double_sig(wilson_matrix_g5_herm(wm2), q.RngState("seed-sig-mul"))
q.displayln_info(f"CHECK: wilson_matrix_g5_herm: {v1:.10f} {v2:.10f}")

q.check_all_files_crc32_info("results")

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
