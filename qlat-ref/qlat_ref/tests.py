"""
Runnable self-checks of the :mod:`qlat_ref.gauge` module.\n
Run with::\n
    cd qlat-ref && python3 -m qlat_ref.tests\n
The checks are deliberately small and readable; they use tiny lattices so
that the pure-NumPy reference implementation runs in a few seconds.
"""

import numpy as np

from qlat_ref import gauge as g

TOTAL_SITE = (2, 2, 2, 2)
SEED = 20240517
EPS = 1e-10
ALPHA = 0.1

def get_random_field(seed=SEED, total_site=TOTAL_SITE):
    return g.get_random_gauge_field(total_site, seed=seed, scale=1.0)

def fix_and_check(fix, U, **kwargs):
    """Fix ``U`` and assert the post-conditions shared by every successful fixing."""
    kwargs.setdefault("max_iter", 10000)
    kwargs.setdefault("eps", EPS)
    U_fixed, info = fix(U, **kwargs)
    assert info["converged"]
    assert info["divergence"] <= EPS
    assert g.check_unitarity(U_fixed, eps=1e-10)
    return U_fixed, info

def check_random_field_is_unitary():
    U = get_random_field()
    assert g.check_unitarity(U, eps=1e-12)
    print(f"check_random_field_is_unitary: max defect = {g.get_unitarity_error(U):.3e}")

def check_trivially_fixed_fields():
    # The unit field is gauge fixed by construction.
    U = g.get_unit_gauge_field(TOTAL_SITE)
    assert g.get_gauge_functional(U) == 1.0
    assert g.get_gauge_divergence(U) < 1e-14
    U_fixed, _ = fix_and_check(g.fix_landau_gauge, U, max_iter=10)
    assert np.max(np.abs(U_fixed - U)) == 0.0
    print("check_unit_field_is_already_fixed: functional = 1, divergence = 0")
    # On a 1^4 lattice a gauge transformation acts at a single site only and
    # cancels in Re Tr[g U g^dagger], so every field is trivially "fixed" and
    # the local gradient must vanish.
    U = g.get_random_gauge_field((1, 1, 1, 1), seed=SEED, scale=1.0)
    functional_before = g.get_gauge_functional(U)
    divergence_before = g.get_gauge_divergence(U)
    assert divergence_before < 1e-14
    U_fixed, info = fix_and_check(g.fix_landau_gauge, U, max_iter=10)
    assert info["n_iter"] == 1
    assert abs(g.get_gauge_functional(U_fixed) - functional_before) < 1e-14
    print(
        "check_single_site_field_is_already_fixed:"
        f" divergence = {divergence_before:.3e}"
    )

def check_landau_gauge_fixing():
    U = get_random_field()
    functional_before = g.get_gauge_functional(U)
    divergence_before = g.get_gauge_divergence(U)
    assert divergence_before > 1e-2
    U_fixed, info = fix_and_check(g.fix_landau_gauge, U, alpha=ALPHA)
    divergence_after = g.get_gauge_divergence(U_fixed)
    assert g.get_gauge_functional(U_fixed) > functional_before
    # The functional must not decrease during the plain steepest-descent
    # iteration (up to the rounding noise of the Jacobi sweep).
    increments = np.diff(info["functional_list"])
    assert np.min(increments) >= -1e-12
    print(
        "check_landau_gauge_fixing:"
        f" functional {functional_before:.10f} -> {info['functional']:.10f},"
        f" divergence {divergence_before:.3e} -> {divergence_after:.3e},"
        f" {info['n_iter']} sweeps"
    )

def check_coulomb_gauge_fixing():
    U = get_random_field(seed=SEED + 1)
    functional_before = g.get_coulomb_functional_list(U)
    U_fixed, info = fix_and_check(g.fix_coulomb_gauge, U, alpha=ALPHA)
    functional_after = g.get_coulomb_functional_list(U_fixed)
    assert all(
        after > before for before, after in zip(functional_before, functional_after)
    )
    # Every time slice is an independent 3D Landau fixing and must satisfy the
    # Coulomb gauge condition on its own.  The convergence criterion uses the
    # RMS divergence over the whole lattice, so a single slice is allowed a
    # factor sqrt(Nt) more.
    divergence_list = g.get_coulomb_divergence_list(U_fixed)
    assert max(divergence_list) <= EPS * np.sqrt(len(divergence_list))
    print(
        "check_coulomb_gauge_fixing:"
        f" {len(functional_after)} time slices,"
        f" functional {np.mean(functional_before):.10f} -> {info['functional']:.10f},"
        f" {info['n_iter']} sweeps"
    )

def check_gauge_covariance():
    # A global gauge transformation g0 (the same SU(3) matrix on every site)
    # leaves the functional Re Tr U_mu(x) invariant, and the plain local
    # iteration is exactly covariant under it: fixing g0 U g0^dagger must give
    # g0 U_fixed g0^dagger and the same functional value.
    U = get_random_field(seed=SEED + 2)
    g_transform = g.get_random_gauge_transform((1, 1, 1, 1), seed=SEED + 3, scale=0.5)
    g_global = np.broadcast_to(g_transform, TOTAL_SITE + (g.NC, g.NC)).copy()
    U_transformed = g.apply_gauge_transform(U, g_global)
    U_fixed, info = fix_and_check(g.fix_landau_gauge, U)
    U_fixed_transformed, info_transformed = fix_and_check(
        g.fix_landau_gauge, U_transformed
    )
    expected = g.apply_gauge_transform(U_fixed, g_global)
    max_difference = float(np.max(np.abs(U_fixed_transformed - expected)))
    assert max_difference < 1e-10
    functional_difference = abs(info["functional"] - info_transformed["functional"])
    assert functional_difference < 1e-12
    print(
        "check_gauge_covariance:"
        f" max |U_fixed^g0 - g0 U_fixed g0^dagger| = {max_difference:.3e},"
        f" |dF| = {functional_difference:.3e}"
    )

def check_over_relaxation_agrees():
    # The converged functional must not depend on the step size (only the
    # number of sweeps does).
    U = get_random_field(seed=SEED + 4)
    U_fixed_a, info_a = fix_and_check(g.fix_landau_gauge, U, alpha=0.05)
    U_fixed_b, info_b = fix_and_check(g.fix_landau_gauge, U, alpha=0.1)
    difference = abs(info_a["functional"] - info_b["functional"])
    assert difference < 1e-10
    print(
        "check_over_relaxation_agrees:"
        f" |F(alpha=0.05) - F(alpha=0.1)| = {difference:.3e}"
    )

def main():
    check_random_field_is_unitary()
    check_trivially_fixed_fields()
    check_landau_gauge_fixing()
    check_coulomb_gauge_fixing()
    check_gauge_covariance()
    check_over_relaxation_agrees()
    print("qlat_ref self-checks: all passed.")

if __name__ == "__main__":
    main()
