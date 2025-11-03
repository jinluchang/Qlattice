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

f_in = q.FieldComplexD(geo, 4 * ls)
f_in.set_rand_g(rs.split("f_in"), 0.0, 1.0)

q.json_results_append("f_in data sig", q.get_data_sig_arr(f_in, rs, 3), check_eps)

q.displayln_info(f"CHECK: ", gf1_qed.geo.expansion_left)
q.displayln_info(f"CHECK: ", gf1_qed.geo.expansion_right)
q.displayln_info(f"CHECK: ", q.geo_resize(geo, 1, 0).show())

f_out = q.multiply_m_dwf_qed(f_in, gf1_qed, mass, m5, ls, is_dagger=False)
q.json_results_append("f_out data sig", q.get_data_sig_arr(f_out, rs, 3), check_eps)

f_out = q.multiply_m_dwf_qed(f_in, gf1_qed, mass, m5, ls, is_dagger=True)
q.json_results_append("f_out dagger data sig", q.get_data_sig_arr(f_out, rs, 3), check_eps)

q.check_log_json(__file__, check_eps=1e-10)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
