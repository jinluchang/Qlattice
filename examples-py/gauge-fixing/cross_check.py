#!/usr/bin/env python3

# Cross-check of the qlat C++/Cython gauge-fixing implementation in this
# example project against the pure-NumPy reference package `qlat-ref`.
#
# This script is NOT part of the CI test suite: the `qlat-ref` sources are a
# sibling of `examples-py/` in the repository, while the nix example build
# only packages `examples-py/`.  Run it manually from a checkout with the nix
# environment and the built extension:
#
#   source result-py-local/bin/setenv-qlat.sh
#   PYTHONPATH=<build install>/lib/python3/dist-packages:$PYTHONPATH \
#     mpiexec -n 2 --oversubscribe --bind-to none \
#     python3 -m mpi4py examples-py/gauge-fixing/cross_check.py
#
# It builds one deterministic random gauge field with `qlat-ref`, fixes it
# with both implementations (same max_iter/eps/alpha) and compares the
# converged functional, the RMS divergence, the unitarity error and (on the
# locally owned sites of each MPI rank) the fixed link matrices.

import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OMP_DYNAMIC"] = "false"

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
# The script directory contains the *sources* of the extension package; drop
# it from sys.path so that the built/installed package is imported instead.
sys.path = [path for path in sys.path if os.path.abspath(path or ".") != SCRIPT_DIR]
sys.path.insert(0, os.path.join(REPO_ROOT, "qlat-ref"))

import qlat as q
import qlat_ref as qr

import qlat_gauge_fixing_ext as qg

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [2, 2, 2, 2],
]

TOTAL_SITE = (2, 2, 2, 2)
SEED = 731
SCALE = 0.7
MAX_ITER = 10000
EPS = 1.0e-10
ALPHA = 0.1

TOLERANCE_FUNCTIONAL = 1.0e-8
TOLERANCE_FIELD = 1.0e-6

def get_buffer_layout(U):
    """Global NumPy field -> (Nt,Nx,Ny,Nz,4,3,3), the qlat buffer layout."""
    return np.transpose(U, (1, 2, 3, 4, 0, 5, 6))

def get_local_buffer(arr, U, geo):
    """Pick the locally owned part of a global NumPy field."""
    xg_arr = geo.xg_arr()
    U_layout = get_buffer_layout(U)
    return U_layout[xg_arr[:, 0], xg_arr[:, 1], xg_arr[:, 2], xg_arr[:, 3]]

def gf_from_numpy(geo, U):
    """Create a qlat GaugeField from a global NumPy field (on every rank)."""
    gf = q.GaugeField(geo)
    arr = np.asarray(gf)
    assert arr.shape == (geo.local_volume, 4, 3, 3)
    arr[...] = get_local_buffer(arr, U, geo)
    return gf

def get_local_difference(gf, U, geo):
    """Max |gf - U| over the locally owned sites."""
    arr = np.asarray(gf)
    return float(np.max(np.abs(arr - get_local_buffer(arr, U, geo))))

def check_gauge(qr_fix, qg_fix, geo, fname):
    U = qr.get_random_gauge_field(TOTAL_SITE, seed=SEED, scale=SCALE)
    gf = gf_from_numpy(geo, U)
    U_ref, ref_info = qr_fix(U, max_iter=MAX_ITER, eps=EPS, alpha=ALPHA)
    gf_info = qg_fix(gf, max_iter=MAX_ITER, eps=EPS, alpha=ALPHA)
    gauge = ref_info["gauge"]
    functional_difference = abs(gf_info["functional"] - ref_info["functional"])
    local_field_difference = get_local_difference(gf, U_ref, geo)
    unitarity_ref = qr.get_unitarity_error(U_ref)
    unitarity_qlat = qg.get_unitarity_error(gf)
    assert ref_info["converged"]
    assert gf_info["converged"]
    assert ref_info["divergence"] <= EPS
    assert gf_info["divergence"] <= EPS
    assert functional_difference <= TOLERANCE_FUNCTIONAL
    assert local_field_difference <= TOLERANCE_FIELD
    assert unitarity_ref <= 1.0e-12
    assert unitarity_qlat <= 1.0e-12
    q.displayln_info(
        f"{fname}: {gauge}: F_ref = {ref_info['functional']:.15e};"
        f" F_qlat = {gf_info['functional']:.15e};"
        f" |dF| = {functional_difference:.3e}"
    )
    q.displayln_info(
        f"{fname}: {gauge}: theta_ref = {ref_info['divergence']:.3e};"
        f" theta_qlat = {gf_info['divergence']:.3e}"
    )
    q.displayln_info(
        f"{fname}: {gauge}: max|U_ref - U_qlat| (local sites) ="
        f" {local_field_difference:.3e}; unitarity_ref = {unitarity_ref:.3e};"
        f" unitarity_qlat = {unitarity_qlat:.3e}"
    )

def main():
    fname = q.get_fname()
    q.begin_with_mpi(size_node_list)
    geo = q.Geometry(q.Coordinate(list(TOTAL_SITE)))
    q.displayln_info(
        f"{fname}: total_site = {TOTAL_SITE}; num_node = {q.get_num_node()};"
        f" local_volume = {geo.local_volume}"
    )
    check_gauge(qr.fix_landau_gauge, qg.fix_landau_gauge, geo, fname)
    check_gauge(qr.fix_coulomb_gauge, qg.fix_coulomb_gauge, geo, fname)
    q.end_with_mpi()
    q.displayln_info("CHECK: cross-check finished successfully.")

if __name__ == "__main__":
    main()
