"""
Module ``qlat.propagator_utils``
=================================\n
Propagator helpers that do not call C++ functions directly: random-U(1)
propagator generation and free-scalar inversion in position space.\n
"""

import qlat_utils as q

from .propagator import (
        mk_rand_u1_src,
        get_rand_u1_sol,
        free_scalar_mom_invert,
        free_scalar_deriv_mom,
        )
from .field_utils_utils import mk_fft

@q.timer_verbose
def mk_rand_u1_prop(inv, sel, rs):
    """
    interface function
    return s_prop
    sel can be psel or fsel
    """
    prop_src, fu1 = mk_rand_u1_src(sel, rs)
    prop_sol = inv * prop_src
    return get_rand_u1_sol(prop_sol, fu1, sel)

@q.timer
def free_scalar_invert(src, mass, *, momtwist=None, mode_fft=1):
    fft_f = mk_fft(is_forward=True, is_normalizing=True, mode_fft=mode_fft)
    fft_b = mk_fft(is_forward=False, is_normalizing=True, mode_fft=mode_fft)
    f = fft_f * src
    free_scalar_mom_invert(f, mass, momtwist)
    sol = fft_b * f
    return sol

@q.timer
def free_scalar_invert_deriv(src, mass, *, momtwist=None, mode_fft=1, deriv=None):
    """
    Free scalar inverse with a lattice derivative, in position space.\n
    Transforms `src` to momentum space, applies the bare derivative factor
    `free_scalar_deriv_mom` with the orders `deriv`, applies the free scalar
    inverse `free_scalar_mom_invert`, and transforms back.  The derivative and
    the inverse commute, so this is the derivative of the free scalar inverse
    (equivalently the free scalar inverse of the derivative source).\n
    `deriv=None` is equivalent to `[0, 0, 0, 0]`, in which case the result
    equals `free_scalar_invert(src, mass, ...)`.
    """
    fft_f = mk_fft(is_forward=True, is_normalizing=True, mode_fft=mode_fft)
    fft_b = mk_fft(is_forward=False, is_normalizing=True, mode_fft=mode_fft)
    f = fft_f * src
    free_scalar_deriv_mom(f, deriv, momtwist)
    free_scalar_mom_invert(f, mass, momtwist)
    sol = fft_b * f
    return sol
