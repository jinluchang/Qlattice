"""
Package ``qlat_ref``
====================\n
Pure-NumPy reference implementations used to cross-check the qlat C++/Cython
code.  The gauge-fixing module is :mod:`qlat_ref.gauge`; its documentation is
``docs/qlat-ref/qlat_gauge.md``.
"""

from .gauge import (
    NC,
    apply_gauge_transform,
    check_unitarity,
    fix_coulomb_gauge,
    fix_landau_gauge,
    get_coulomb_divergence_list,
    get_coulomb_functional_list,
    get_dirs,
    get_divergence_from_gradient,
    get_gauge_divergence,
    get_gauge_functional,
    get_local_gauge_gradient,
    get_matrix_adjoint,
    get_random_gauge_field,
    get_random_gauge_transform,
    get_staple_sum,
    get_su3_exp,
    get_traceless_antihermitian_part,
    get_unit_gauge_field,
    get_unitarity_error,
)

__all__ = [
    "NC",
    "apply_gauge_transform",
    "check_unitarity",
    "fix_coulomb_gauge",
    "fix_landau_gauge",
    "get_coulomb_divergence_list",
    "get_coulomb_functional_list",
    "get_dirs",
    "get_divergence_from_gradient",
    "get_gauge_divergence",
    "get_gauge_functional",
    "get_local_gauge_gradient",
    "get_matrix_adjoint",
    "get_random_gauge_field",
    "get_random_gauge_transform",
    "get_staple_sum",
    "get_su3_exp",
    "get_traceless_antihermitian_part",
    "get_unit_gauge_field",
    "get_unitarity_error",
]
