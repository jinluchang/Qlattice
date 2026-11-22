#!/usr/bin/env python3

# Plain local (no Fourier acceleration) Landau and Coulomb gauge fixing test
# for the qlat C++/Cython gauge-fixing extension of this example project.
#
# Run through the examples-py Makefile (target gauge-fixing-ext.log), which
# builds the meson project and launches this script with 2 MPI processes.

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OMP_DYNAMIC"] = "false"

import qlat as q

import qlat_gauge_fixing_ext as qg

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [2, 2, 2, 2],
]

TOTAL_SITE = [4, 4, 4, 4]
EPS = 1.0e-8
ALPHA = 0.1
MAX_ITER = 5000

def show_info(fname, gf, tag):
    plaq = gf.plaq()
    link_trace = gf.link_trace()
    q.displayln_info(
        f"{fname}: {tag}: plaq = {plaq:.15e}; link_trace = {link_trace:.15e}"
    )

def check_landau_gauge_fixing(fname, geo, rs):
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split("gf-landau"), 0.5, 1)
    show_info(fname, gf, "landau before")
    plaq_before = gf.plaq()
    functional_before = qg.get_gauge_functional(gf, [0, 1, 2, 3])
    divergence_before = qg.get_gauge_divergence(gf, [0, 1, 2, 3])
    info = qg.fix_landau_gauge(gf, max_iter=MAX_ITER, eps=EPS, alpha=ALPHA)
    show_info(fname, gf, "landau after")
    plaq_after = gf.plaq()
    unitarity_error = qg.get_unitarity_error(gf)
    q.displayln_info(
        f"{fname}: landau: functional {functional_before:.15e} ->"
        f" {info['functional']:.15e} in {info['n_iter']} sweeps;"
        f" divergence {divergence_before:.3e} -> {info['divergence']:.3e};"
        f" unitarity error = {unitarity_error:.3e}"
    )
    assert info["gauge"] == "landau"
    assert info["converged"]
    assert info["divergence"] <= EPS
    assert info["functional"] > functional_before
    assert abs(plaq_after - plaq_before) < 1.0e-10
    assert unitarity_error < 1.0e-10
    q.json_results_append(f"{fname}: landau plaq before", plaq_before, 1.0e-10)
    q.json_results_append(f"{fname}: landau plaq after", plaq_after, 1.0e-10)
    q.json_results_append(
        f"{fname}: landau functional before", functional_before, 1.0e-10
    )
    q.json_results_append(f"{fname}: landau functional", info["functional"], 1.0e-10)

def check_coulomb_gauge_fixing(fname, geo, rs):
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split("gf-coulomb"), 0.5, 1)
    show_info(fname, gf, "coulomb before")
    plaq_before = gf.plaq()
    functional_before = qg.get_gauge_functional(gf, [1, 2, 3])
    divergence_before = qg.get_gauge_divergence(gf, [1, 2, 3])
    info = qg.fix_coulomb_gauge(gf, max_iter=MAX_ITER, eps=EPS, alpha=ALPHA)
    show_info(fname, gf, "coulomb after")
    plaq_after = gf.plaq()
    unitarity_error = qg.get_unitarity_error(gf)
    q.displayln_info(
        f"{fname}: coulomb: functional {functional_before:.15e} ->"
        f" {info['functional']:.15e} in {info['n_iter']} sweeps;"
        f" divergence {divergence_before:.3e} -> {info['divergence']:.3e};"
        f" unitarity error = {unitarity_error:.3e}"
    )
    assert info["gauge"] == "coulomb"
    assert info["converged"]
    assert info["divergence"] <= EPS
    assert info["functional"] > functional_before
    assert abs(plaq_after - plaq_before) < 1.0e-10
    assert unitarity_error < 1.0e-10
    q.json_results_append(f"{fname}: coulomb plaq before", plaq_before, 1.0e-10)
    q.json_results_append(f"{fname}: coulomb plaq after", plaq_after, 1.0e-10)
    q.json_results_append(
        f"{fname}: coulomb functional before", functional_before, 1.0e-10
    )
    q.json_results_append(f"{fname}: coulomb functional", info["functional"], 1.0e-10)

def check_distributed_field(fname, geo):
    gf = q.GaugeField(geo)
    rs = q.RngState("gauge-fixing-distributed")
    gf.set_rand(rs.split("gf"), 0.5, 1)
    info = qg.fix_landau_gauge(gf, max_iter=MAX_ITER, eps=EPS, alpha=ALPHA)
    assert info["converged"]
    assert qg.get_gauge_divergence(gf, [0, 1, 2, 3]) <= EPS
    q.displayln_info(
        f"{fname}: distributed: id_node = {q.get_id_node()};"
        f" num_node = {q.get_num_node()}; local_volume = {geo.local_volume};"
        f" total_volume = {geo.total_volume}"
    )
    q.json_results_append(f"{fname}: num_node = {q.get_num_node()}")
    q.json_results_append(f"{fname}: local_volume = {geo.local_volume}")
    q.json_results_append(f"{fname}: total_volume = {geo.total_volume}")

def main():
    fname = q.get_fname()
    q.begin_with_mpi(size_node_list)
    total_site = q.Coordinate(TOTAL_SITE)
    geo = q.Geometry(total_site)
    q.json_results_append(f"{fname}: geo = {geo.show()}")
    rs = q.RngState("gauge-fixing")
    check_landau_gauge_fixing(fname, geo, rs)
    check_coulomb_gauge_fixing(fname, geo, rs)
    check_distributed_field(fname, geo)
    q.timer_display()
    if q.is_test():
        q.check_log_json(__file__, check_eps=1.0e-10)
    q.end_with_mpi()
    q.displayln_info("CHECK: finished successfully.")

if __name__ == "__main__":
    main()
