#!/usr/bin/env python3

"""
Tests for qlat.selected_points_io.convert_selected_points_dist_type.

Covers:
  1. l -> g conversion (dist type changes, data accessible)
  2. g -> l conversion (dist type changes)
  3. l -> g -> l round-trip (compare g-type signatures via intermediate g)
  4. g -> l -> g round-trip (compare g-type signatures)
  5. sp_points_dist_type assertion (correct / incorrect)
  6. root parameter (explicit root)
  7. Invalid conversions raise exceptions
  8. Multiple SelectedPoints types and multiplicities
"""

import numpy as np
import qlat as q

q.begin_with_mpi()

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
multiplicity = 3

rs = q.RngState("seed convert_selected_points_dist_type")

# ---------- helpers ----------

def make_local_psel(geo, rs):
    """Create a PointsSelection with points_dist_type='l' (local) from FieldSelection."""
    fsel = q.FieldSelection(geo)
    fsel.set_rand(total_site, 1, rs.split("fsel"))
    psel = q.PointsSelection(fsel)
    return psel

def make_global_psel(geo, n_points, rs):
    """Create a PointsSelection with points_dist_type='g' (global)."""
    psel = q.PointsSelection()
    psel.set_rand(total_site, n_points, rs.split("psel"))
    return psel

# ---------- test 1: l -> g conversion ----------

def test_l_to_g(geo, multiplicity, seed):
    """Convert l -> g. Verify dist type changes and data is accessible."""
    fname = q.get_fname()
    rs = q.RngState(f"seed {fname} {seed}")
    q.json_results_append(f"{fname}: total_site={total_site}, multiplicity={multiplicity}")

    psel_l = make_local_psel(geo, rs.split("psel_l"))
    q.json_results_append(f"psel_l.points_dist_type={psel_l.points_dist_type}")

    sp_l = q.SelectedPointsComplexD(psel_l, multiplicity)
    sp_l.set_rand(rs.split("sp_l"))
    q.json_results_append(f"sp_l.points_dist_type={sp_l.points_dist_type}")

    # l -> g (root=None -> root=0 with broadcast)
    sp_g = q.convert_selected_points_dist_type(sp_l, "g")
    q.json_results_append(f"sp_g.points_dist_type={sp_g.points_dist_type}")

    assert sp_g.points_dist_type == "g"
    # After l->g with default root, all nodes have the full data
    data = np.asarray(sp_g)
    assert data.size > 0, "l->g data is empty"
    q.json_results_append("l->g conversion: PASS")

# ---------- test 2: g -> l conversion ----------

def test_g_to_l(geo, multiplicity, seed):
    """Convert g -> l. Verify dist type changes."""
    fname = q.get_fname()
    rs = q.RngState(f"seed {fname} {seed}")
    q.json_results_append(f"{fname}: total_site={total_site}, multiplicity={multiplicity}")

    n_points = total_site.volume() // 16
    psel_g = make_global_psel(geo, n_points, rs.split("psel_g"))
    q.json_results_append(f"psel_g.points_dist_type={psel_g.points_dist_type}")

    sp_g = q.SelectedPointsComplexD(psel_g, multiplicity)
    sp_g.set_rand(rs.split("sp_g"))
    q.json_results_append(f"sp_g.points_dist_type={sp_g.points_dist_type}")

    # g -> l
    sp_l = q.convert_selected_points_dist_type(sp_g, "l")
    q.json_results_append(f"sp_l.points_dist_type={sp_l.points_dist_type}")

    assert sp_l.points_dist_type == "l"
    q.json_results_append("g->l conversion: PASS")

# ---------- test 3: l -> g -> l round-trip (g sig comparison) ----------

def test_l_to_g_roundtrip(geo, multiplicity, seed):
    """
    l -> g -> l round-trip.
    Compare the g-type signatures before and after the round-trip.
    We go through g for both comparisons to avoid psel-ordering issues.
    """
    fname = q.get_fname()
    rs = q.RngState(f"seed {fname} {seed}")
    q.json_results_append(f"{fname}: total_site={total_site}, multiplicity={multiplicity}")

    psel_l = make_local_psel(geo, rs.split("psel_l"))

    # Create a g-type psel with the same points for signature comparison
    psel_g = make_global_psel(geo, psel_l.n_points, rs.split("psel_g_ref"))
    sp_g_ref = q.SelectedPointsComplexD(psel_g, multiplicity)
    sp_g_ref.set_rand(rs.split("sp_g_ref"))
    sig_g_ref = q.get_data_sig_arr(sp_g_ref, rs.split("sig"), 2)

    # Now create local sp from l->g conversion of the global sp
    sp_l = q.convert_selected_points_dist_type(sp_g_ref, "l")
    q.json_results_append(f"sp_l.points_dist_type={sp_l.points_dist_type}")

    # l -> g round-trip
    sp_g2 = q.convert_selected_points_dist_type(sp_l, "g")
    q.json_results_append(f"sp_g2.points_dist_type={sp_g2.points_dist_type}")

    sig_g2 = q.get_data_sig_arr(sp_g2, rs.split("sig"), 2)

    # Both are "g" type with the same rng, so signatures should match
    assert np.allclose(sig_g_ref, sig_g2, atol=1e-10), (
        f"l->g->l sig mismatch: {sig_g_ref} vs {sig_g2}"
    )
    q.json_results_append("l->g->l round-trip: PASS")

# ---------- test 4: g -> l -> g round-trip ----------

def test_g_to_l_roundtrip(geo, multiplicity, seed):
    """g -> l -> g round-trip. Compare g-type signatures."""
    fname = q.get_fname()
    rs = q.RngState(f"seed {fname} {seed}")
    q.json_results_append(f"{fname}: total_site={total_site}, multiplicity={multiplicity}")

    n_points = total_site.volume() // 16
    psel_g = make_global_psel(geo, n_points, rs.split("psel_g"))
    sp_g = q.SelectedPointsComplexD(psel_g, multiplicity)
    sp_g.set_rand(rs.split("sp_g"))

    sig_g = q.get_data_sig_arr(sp_g, rs.split("sig"), 2)
    q.json_results_append("g sig", sig_g, 1e-10)

    # g -> l -> g
    sp_l = q.convert_selected_points_dist_type(sp_g, "l")
    sp_g2 = q.convert_selected_points_dist_type(sp_l, "g")

    sig_g2 = q.get_data_sig_arr(sp_g2, rs.split("sig"), 2)
    q.json_results_append("g2 sig", sig_g2, 1e-10)

    assert np.allclose(sig_g, sig_g2, atol=1e-10), (
        f"g->l->g sig mismatch: {sig_g} vs {sig_g2}"
    )
    q.json_results_append("g->l->g round-trip: PASS")

# ---------- test 5: sp_points_dist_type assertion ----------

def test_sp_points_dist_type_assertion(geo, multiplicity, seed):
    fname = q.get_fname()
    rs = q.RngState(f"seed {fname} {seed}")
    q.json_results_append(f"{fname}: total_site={total_site}")

    n_points = total_site.volume() // 16
    psel = make_global_psel(geo, n_points, rs.split("psel"))
    sp = q.SelectedPointsComplexD(psel, multiplicity)
    sp.set_rand(rs.split("sp"))

    # Should pass when correct type specified
    sp_g = q.convert_selected_points_dist_type(sp, "l", sp_points_dist_type="g")
    assert sp_g.points_dist_type == "l"
    q.json_results_append("sp_points_dist_type assertion (correct): PASS")

    # Should fail when wrong type specified
    try:
        q.convert_selected_points_dist_type(sp, "l", sp_points_dist_type="l")
        q.json_results_append(
            "sp_points_dist_type assertion (wrong): FAIL - no exception raised"
        )
    except AssertionError:
        q.json_results_append(
            "sp_points_dist_type assertion (wrong): PASS (raised AssertionError)"
        )

# ---------- test 6: root parameter ----------

def test_root_parameter(geo, multiplicity, seed):
    fname = q.get_fname()
    rs = q.RngState(f"seed {fname} {seed}")
    q.json_results_append(f"{fname}: total_site={total_site}")

    psel_l = make_local_psel(geo, rs.split("psel"))
    sp_l = q.SelectedPointsComplexD(psel_l, multiplicity)
    sp_l.set_rand(rs.split("sp"))

    # l -> g with explicit root=0 (no broadcast)
    sp_g = q.convert_selected_points_dist_type(sp_l, "g", root=0)
    assert sp_g.points_dist_type == "g"

    # g -> l with explicit root=0
    sp_l2 = q.convert_selected_points_dist_type(sp_g, "l", root=0)
    assert sp_l2.points_dist_type == "l"

    q.json_results_append("root parameter round-trip: PASS")

# ---------- test 7: invalid conversions raise ----------

def test_invalid_conversion(geo, multiplicity, seed):
    fname = q.get_fname()
    rs = q.RngState(f"seed {fname} {seed}")
    q.json_results_append(f"{fname}: total_site={total_site}")

    n_points = total_site.volume() // 16

    # l -> l should raise
    psel_l = make_local_psel(geo, rs.split("psel"))
    sp_l = q.SelectedPointsComplexD(psel_l, multiplicity)
    try:
        q.convert_selected_points_dist_type(sp_l, "l")
        q.json_results_append("l->l invalid: FAIL (no exception)")
    except Exception:
        q.json_results_append("l->l invalid: PASS (raised exception)")

    # g -> g should raise
    psel_g = make_global_psel(geo, n_points, rs.split("psel_g"))
    sp_g = q.SelectedPointsComplexD(psel_g, multiplicity)
    try:
        q.convert_selected_points_dist_type(sp_g, "g")
        q.json_results_append("g->g invalid: FAIL (no exception)")
    except Exception:
        q.json_results_append("g->g invalid: PASS (raised exception)")

# ---------- test 8: multiple types and multiplicities ----------

def test_types_and_multiplicity(geo, seed):
    fname = q.get_fname()
    rs = q.RngState(f"seed {fname} {seed}")
    q.json_results_append(f"{fname}: total_site={total_site}")

    n_points = total_site.volume() // 16

    for cls_name, cls in [
        ("SelectedPointsRealD", q.SelectedPointsRealD),
        ("SelectedPointsRealF", q.SelectedPointsRealF),
        ("SelectedPointsComplexD", q.SelectedPointsComplexD),
        ("SelectedPointsComplexF", q.SelectedPointsComplexF),
    ]:
        for mult in [1, 2, 5]:
            # g -> l -> g round-trip
            psel = make_global_psel(geo, n_points, rs.split(f"psel_{cls_name}_{mult}"))
            sp = cls(psel, mult)
            sp.set_rand(rs.split(f"sp_{cls_name}_{mult}"))

            sig_before = q.get_data_sig_arr(
                sp, rs.split(f"sig_{cls_name}_{mult}"), 2
            )

            sp_l = q.convert_selected_points_dist_type(sp, "l")
            assert sp_l.points_dist_type == "l"
            sp_g2 = q.convert_selected_points_dist_type(sp_l, "g")
            assert sp_g2.points_dist_type == "g"

            sig_after = q.get_data_sig_arr(
                sp_g2, rs.split(f"sig_{cls_name}_{mult}"), 2
            )
            assert np.allclose(sig_before, sig_after, atol=1e-12), (
                f"{cls_name}(mult={mult}): g->l->g signature mismatch"
            )
            q.json_results_append(f"{cls_name}(mult={mult}): PASS")

# ---------- run all tests ----------

for seed in range(1):
    test_l_to_g(geo, multiplicity, seed)
    test_g_to_l(geo, multiplicity, seed)
    test_l_to_g_roundtrip(geo, multiplicity, seed)
    test_g_to_l_roundtrip(geo, multiplicity, seed)
    test_types_and_multiplicity(geo, seed)
    test_sp_points_dist_type_assertion(geo, multiplicity, seed)
    test_root_parameter(geo, multiplicity, seed)
    test_invalid_conversion(geo, multiplicity, seed)

q.json_results_append("convert_selected_points_dist_type: all tests passed")
q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
