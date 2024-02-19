#!/usr/bin/env python3

json_results = []

import qlat as q
import numpy as np

q.begin_with_mpi()

has_cuba = q.has_cuba()

q.displayln_info(f"CHECK: has_cuba={has_cuba}")

if not has_cuba:
    q.displayln_info(f"CHECK: quit due to not having CUBA.")
    q.end_with_mpi()
    exit()

q.test_integration_multi_dimensional()

cx = q.Coordinate([ 1, 2, 3, 4, ])
cy = q.Coordinate([ 3, 1, 2, 1, ])
cz = q.Coordinate([ 0, 0, 0, 0, ])
total_site = q.Coordinate([ 16, 16, 16, 32, ])
muon_mass = 0.4

c1 = q.CoordinateD([ 1, 2, 3, 4, ]) * 0.4
c2 = q.CoordinateD([ 3, 1, 2, 1, ]) * 0.4
c3 = q.CoordinateD([ 0, 0, 0, 0, ]) * 0.4

q.displayln_info(0, f"{c1, c2, c3}")

epsabs = 1e-8
epsrel = 1e-3

ve = None

if True:
    fname = "calc_muon_line_m"
    v0 = q.calc_muon_line_m(c1, c2, epsabs * 10, epsrel * 10)
    json_results.append((f"{fname}: v0 sig", q.get_double_sig(v0, q.RngState()),))
    v1 = q.calc_muon_line_m(c1, c2, epsabs, epsrel)
    json_results.append((f"{fname}: v1 sig", q.get_double_sig(v1, q.RngState()),))
    ve = v1
    q.displayln_info(np.sqrt(q.qnorm(v0 - ve) / q.qnorm(ve)))
    if False:
        v2 = q.calc_muon_line_m(c1, c2, epsabs / 10, epsrel / 10)
        json_results.append((f"{fname}: v2 sig", q.get_double_sig(v2, q.RngState()),))
        ve = v2
        q.displayln_info(np.sqrt(q.qnorm(v0 - ve) / q.qnorm(ve)))
        q.displayln_info(np.sqrt(q.qnorm(v1 - ve) / q.qnorm(ve)))

if True:
    q.compute_save_muonline_interpolation(
            f"results/muon-line-interpolation-data/{0:010d}",
            [ 2, 2, 2, 2, 2, ],
            epsabs, epsrel)
    q.load_multiple_muonline_interpolations("results/muon-line-interpolation-data/", [ 0, ])
    fname = "get_muon_line_m"
    vv1 = q.get_muon_line_m(c1, c2, c3, 0, epsabs, epsrel)
    json_results.append((f"{fname}: vv1 sig", q.get_double_sig(vv1, q.RngState()),))
    vv2 = q.get_muon_line_m(c2, c3, c1, 0, epsabs, epsrel)
    json_results.append((f"{fname}: vv2 sig", q.get_double_sig(vv2, q.RngState()),))
    vv3 = q.get_muon_line_m(c3, c1, c2, 0, epsabs, epsrel)
    json_results.append((f"{fname}: vv3 sig", q.get_double_sig(vv3, q.RngState()),))
    vv4 = q.get_muon_line_m(c3, c2, c1, 0, epsabs, epsrel)
    json_results.append((f"{fname}: vv4 sig", q.get_double_sig(vv4, q.RngState()),))
    vv5 = q.get_muon_line_m(c2, c1, c3, 0, epsabs, epsrel)
    json_results.append((f"{fname}: vv5 sig", q.get_double_sig(vv5, q.RngState()),))
    vv6 = q.get_muon_line_m(c1, c3, c2, 0, epsabs, epsrel)
    json_results.append((f"{fname}: vv6 sig", q.get_double_sig(vv6, q.RngState()),))
    q.displayln_info(np.sqrt(q.qnorm(vv1)))
    q.displayln_info(np.sqrt(q.qnorm(vv1 - ve)))
    q.clear_muon_line_interpolations()

if 0 == q.get_id_node():
    import os
    json_fn_name = os.path.splitext(__file__)[0] + ".log.json"
    q.qtouch(json_fn_name + ".new", q.json_dumps(json_results, indent=1))
    if q.does_file_exist_qar(json_fn_name):
        json_results_load = q.json_loads(q.qcat(json_fn_name))
        for i, p in enumerate(json_results_load):
            nl, vl = p
            n, v = json_results[i]
            if n != nl:
                q.displayln(f"CHECK: {i} {p}")
                q.displayln("CHECK: ERROR: JSON results item does not match.")
                assert False
            if abs(v - vl) > 1e-5 * (abs(v) + abs(vl)):
                q.displayln(f"CHECK: {i} {p}")
                q.displayln("CHECK: ERROR: JSON results value does not match.")
        if len(json_results) != len(json_results_load):
            q.displayln(f"CHECK: len(json_results)={len(json_results)} load:{len(json_results_load)}")
            q.displayln("CHECK: ERROR: JSON results len does not match.")

q.timer_display()

q.end_with_mpi()

q.displayln_info("CHECK: finished successfully.")
