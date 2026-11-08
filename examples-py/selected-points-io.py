#!/usr/bin/env python3

"""
Test for ``qlat.selected_points_io.load_selected_points_list``.\n
Distributed read: each file is read by one MPI rank (no duplication).
Shuffle: l_from_g plan redistributes to Local distribution.\n
Verifies:
  1. Single-file: save → load → shuffle to local → check qnorm preservation.
  2. Batch (multi-file): same verification for multiple files with varying sizes.\n
Every logged result is a global aggregate (``q.glb_sum``), so the reference is
independent of the MPI layout; ``verify_local_sp`` also asserts that each node
holds exactly the points that it owns locally. ``size_node_list`` only lists
layouts that divide both 4⁴ and 6⁴, so the test runs on 1, 2, 4, 8 or 16 ranks.
"""

import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 2, 2, 1],
    [1, 2, 2, 2],
    [2, 2, 2, 2],
]

def verify_local_sp(geo, sp_orig, sp_local):
    """
    Verify that a locally-distributed sp preserves the global point count and
    qnorm, and that each node holds only its locally-owned points.
    The l_from_g shuffle sorts points by gindex, so the ordering changes,
    but the total norm (sum of squared values across all points) is invariant.
    """
    n_points = q.glb_sum(sp_local.n_points)
    assert n_points == sp_orig.n_points, (
        f"n_points mismatch: n_points={n_points} orig_n_points={sp_orig.n_points}"
    )
    local_norm_sq = sp_local.qnorm()
    total_norm_sq = q.glb_sum(local_norm_sq)
    orig_norm_sq = sp_orig.qnorm()
    assert abs(total_norm_sq - orig_norm_sq) <= 1e-12 * max(1.0, abs(orig_norm_sq)), (
        f"qnorm mismatch: total_norm_sq={total_norm_sq} orig_norm_sq={orig_norm_sq}"
    )
    if sp_local.n_points > 0:
        for xg in sp_local.psel:
            assert geo.is_local_xg(xg), (
                f"point {xg} on node {q.get_id_node()} is not locally owned"
            )

def test_single(total_site, multiplicity, seed):
    fname = q.get_fname()
    q.json_results_append(f"{fname}: {total_site} {multiplicity} {seed}")
    rs = q.RngState(f"seed {fname} {seed}")
    n_points = total_site.volume() // 16
    geo = q.Geometry(total_site)
    # Create random psel and sp
    psel = q.PointsSelection()
    psel.set_rand(total_site, n_points, rs.split("psel"))
    sp = q.SelectedPointsRealD(psel, multiplicity)
    sp.set_rand(rs.split("sp"), 1.0, -1.0)
    q.json_results_append("sp.qnorm()", sp.qnorm())
    # Save
    path_psel = f"results/test-single-{multiplicity}-{seed}-psel.lati"
    path_sp = f"results/test-single-{multiplicity}-{seed}-sp.lat"
    psel.save(path_psel)
    sp.save(path_sp)
    # Load psel on all nodes
    psel_load = q.PointsSelection()
    psel_load.load(path_psel, is_sync_node=False)
    # Load and shuffle to local
    sp_list = q.load_selected_points_list(q.SelectedPointsRealD, [psel_load], [path_sp])
    assert len(sp_list) == 1
    sp_local = sp_list[0]
    q.json_results_append(f"sp_local.n_points={q.glb_sum(sp_local.n_points)}")
    q.json_results_append("sp_local.qnorm()", q.glb_sum(sp_local.qnorm()))
    # Verify correctness
    verify_local_sp(geo, sp, sp_local)

def test_batch(total_site, multiplicity, seed):
    fname = q.get_fname()
    q.json_results_append(f"{fname}: {total_site} {multiplicity} {seed}")
    rs = q.RngState(f"seed {fname} {seed}")
    geo = q.Geometry(total_site)
    n_files = 3
    psel_list = []
    sp_list_orig = []
    path_psel_list = []
    path_sp_list = []
    for k in range(n_files):
        n_points = total_site.volume() // (16 * (k + 1))
        psel = q.PointsSelection()
        psel.set_rand(total_site, n_points, rs.split(f"psel-{k}"))
        sp = q.SelectedPointsRealD(psel, multiplicity)
        sp.set_rand(rs.split(f"sp-{k}"), 1.0, -1.0)
        path_psel = f"results/test-batch-{multiplicity}-{seed}-psel-{k}.lati"
        path_sp = f"results/test-batch-{multiplicity}-{seed}-sp-{k}.lat"
        psel.save(path_psel)
        sp.save(path_sp)
        psel_list.append(psel)
        sp_list_orig.append(sp)
        path_psel_list.append(path_psel)
        path_sp_list.append(path_sp)
    # Load psels on all nodes
    psel_load_list = []
    for path_psel in path_psel_list:
        psel_load = q.PointsSelection()
        psel_load.load(path_psel, is_sync_node=False)
        psel_load_list.append(psel_load)
    # Load and shuffle to local
    sp_local_list = q.load_selected_points_list(
        q.SelectedPointsRealD, psel_load_list, path_sp_list
    )
    assert len(sp_local_list) == n_files
    # Verify each
    for k in range(n_files):
        sp_local = sp_local_list[k]
        sp_orig = sp_list_orig[k]
        q.json_results_append(f"sp_local[{k}].n_points={q.glb_sum(sp_local.n_points)}")
        q.json_results_append(f"sp_local[{k}].qnorm()", q.glb_sum(sp_local.qnorm()))
        verify_local_sp(geo, sp_orig, sp_local)

# --- test harness below (must be AFTER function definitions) ---

q.begin_with_mpi(size_node_list)

total_site_list = [
    q.Coordinate([4, 4, 4, 4]),
    q.Coordinate([6, 6, 6, 6]),
]
multiplicity_list = [1, 3]

for total_site in total_site_list:
    for multiplicity in multiplicity_list:
        for seed in range(1):
            test_single(total_site, multiplicity, seed)
            test_batch(total_site, multiplicity, seed)

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
