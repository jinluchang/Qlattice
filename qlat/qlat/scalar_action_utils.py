"""
Module ``qlat.scalar_action_utils``
====================================\n
cqlat-compatible entry points for ``ScalarAction`` that do not call C++
functions directly.\n
"""

def hmc_estimate_mass_scalar_action(sa, masses, field_ft, force_ft, phi0):
    """
    cqlat-compatible entry point for ``ScalarAction.hmc_estimate_mass``.
    """
    sa.hmc_estimate_mass(masses, field_ft, force_ft, phi0)

def to_mass_factor_scalar_action(sa, sin_domega):
    """
    cqlat-compatible entry point for ``ScalarAction.to_mass_factor``.
    """
    sa.to_mass_factor(sin_domega)
