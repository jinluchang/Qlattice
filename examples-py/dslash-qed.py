#!/usr/bin/env python3

import qlat as q
import numpy as np

q.begin_with_mpi()

check_eps = 1e-10

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)

rs = q.RngState("seed")

gf_qed_phase = q.FieldRealD(geo, 4)
gf_qed_phase.set_rand_g(rs.split("gf_qed_phase"), 0.0, 1.0)

gf_qed = q.FieldComplexD(geo, 4)
gf_qed[:] = np.exp(1j * gf_qed_phase[:])

gf1_qed = q.mk_left_expanded_field(gf_qed)

mass = 0.1
ls = 16
m5 = 1.0

q.json_results_append("test 5d multiplication and inversion")

f_in = q.FieldComplexD(geo, 4 * ls)
f_in.set_rand_g(rs.split("f_in"), 0.0, 1.0)
q.json_results_append("f_in data sig", q.get_data_sig_arr(f_in, rs, 3), check_eps)

q.displayln_info(f"CHECK: ", gf1_qed.geo.expansion_left)
q.displayln_info(f"CHECK: ", gf1_qed.geo.expansion_right)
q.displayln_info(f"CHECK: ", q.geo_resize(geo, 1, 0).show())

f_out = q.multiply_m_dwf_qed(f_in, gf1_qed, mass, m5, ls, is_dagger=False)
q.json_results_append("f_out data sig", q.get_data_sig_arr(f_out, rs, 3), check_eps)

f_out = q.multiply_m_dwf_qed(f_out, gf1_qed, mass, m5, ls, is_dagger=True)
f_in2 = q.cg_with_m_dwf_qed(f_out, gf1_qed, mass, m5, ls, is_dagger=False, max_num_iter=100)
q.json_results_append("f_in2 data sig", q.get_data_sig_arr(f_in2, rs, 3), 1e-6)

f_in2_diff = f_in2.copy()
f_in2_diff -= f_in
assert f_in2_diff.qnorm() / f_in.qnorm() < 1e-6

f_out = q.multiply_m_dwf_qed(f_in, gf1_qed, mass, m5, ls, is_dagger=True)
q.json_results_append("f_out dagger data sig", q.get_data_sig_arr(f_out, rs, 3), check_eps)

f_out = q.multiply_m_dwf_qed(f_out, gf1_qed, mass, m5, ls, is_dagger=False)
f_in3 = q.cg_with_m_dwf_qed(f_out, gf1_qed, mass, m5, ls, is_dagger=True, max_num_iter=100)
q.json_results_append("f_in3 dagger data sig", q.get_data_sig_arr(f_in3, rs, 3), 1e-6)

f_in3_diff = f_in3.copy()
f_in3_diff -= f_in
assert f_in3_diff.qnorm() / f_in.qnorm() < 1e-6

q.json_results_append("test 4d inversion")

f_in4d = q.FieldComplexD(geo, 4)
f_in4d.set_rand_g(rs.split("f_in4d"), 0.0, 1.0)
q.json_results_append("f_in4d data sig", q.get_data_sig_arr(f_in4d, rs, 3), check_eps)

f_out4d = q.invert_dwf_qed(f_in4d, gf1_qed, mass, m5, ls, is_dagger=False)
q.json_results_append("f_out4d data sig", q.get_data_sig_arr(f_out4d, rs, 3), 1e-6)

f_out4d = q.invert_dwf_qed(f_in4d, gf1_qed, mass, m5, ls, is_dagger=True)
q.json_results_append("f_out4d dagger data sig", q.get_data_sig_arr(f_out4d, rs, 3), 1e-6)


q.json_results_append("test 4d propagator")

sp_in = q.SpinProp(geo, 1)
sp_in.set_rand_g(rs.split("sp_in"), 0.0, 1.0)
q.json_results_append("sp_in data sig", q.get_data_sig_arr(sp_in, rs, 3), check_eps)

sp_out = q.invert_qed(sp_in, gf1_qed, mass, m5, ls, is_dagger=False)
q.json_results_append("sp_out data sig", q.get_data_sig_arr(sp_out, rs, 3), 1e-6)

sp_out = q.invert_qed(sp_in, gf1_qed, mass, m5, ls, is_dagger=True)
q.json_results_append("sp_out dagger data sig", q.get_data_sig_arr(sp_out, rs, 3), 1e-6)

q.json_results_append("test 4d propagator with free gauge field")

gf_qed = q.FieldComplexD(geo, 4)
gf_qed.set_unit()
gf1_qed = q.mk_left_expanded_field(gf_qed)

mass = 0.1
ls = 64
m5 = 1.0

sp_in = q.SpinProp(geo, 1)
sp_in.set_rand_g(rs.split("sp_in"), 0.0, 1.0)
q.json_results_append("sp_in data sig", q.get_data_sig_arr(sp_in, rs, 3), check_eps)

sp_out = q.invert_qed(sp_in, gf1_qed, mass, m5, ls, is_dagger=False, max_num_iter=200)
q.json_results_append("sp_out data sig", q.get_data_sig_arr(sp_out, rs, 3), 1e-6)

sp_out2 = q.free_invert(sp_in, mass, m5)
q.json_results_append("sp_out2 data sig", q.get_data_sig_arr(sp_out2, rs, 3), 1e-6)

sp_out2_diff = sp_out2.copy()
sp_out2_diff -= sp_out
assert sp_out2_diff.qnorm() / sp_out2.qnorm() < 1e-6

q.check_log_json(__file__, check_eps=1e-10)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
