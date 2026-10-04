"""
Module ``qlat.propagator_utils``
=================================\n
Propagator helpers that do not call C++ functions directly: random-U(1)
propagator generation and free-scalar inversion in position space.\n
"""

class q:
    from qlat_utils import (
        timer_verbose,
        timer,
    )
    from .propagator import (
        mk_rand_u1_src,
        get_rand_u1_sol,
        free_scalar_mom_invert,
        free_scalar_deriv_mom,
    )
    from .field_utils_utils import (
        mk_fft,
    )

@q.timer_verbose
def mk_rand_u1_prop(inv, sel, rs):
    """
    interface function
    return s_prop
    sel can be psel or fsel
    """
    prop_src, fu1 = q.mk_rand_u1_src(sel, rs)
    prop_sol = inv * prop_src
    return q.get_rand_u1_sol(prop_sol, fu1, sel)

@q.timer
def free_scalar_invert(src, mass, *, momtwist=None, mode_fft=1):
    fft_f = q.mk_fft(is_forward=True, is_normalizing=True, mode_fft=mode_fft)
    fft_b = q.mk_fft(is_forward=False, is_normalizing=True, mode_fft=mode_fft)
    f = fft_f * src
    q.free_scalar_mom_invert(f, mass, momtwist)
    sol = fft_b * f
    return sol

@q.timer
def free_scalar_invert_deriv(
    src, mass, *, momtwist=None, mode_fft=1, deriv=None, even_deriv_kernel=None
):
    """
    Free scalar inverse with a lattice derivative, in position space.\n
    Transforms `src` to momentum space, applies the bare derivative factor
    `free_scalar_deriv_mom` with the orders `deriv`, applies the free scalar
    inverse `free_scalar_mom_invert`, and transforms back.  The derivative and
    the inverse commute, so this is the derivative of the free scalar inverse
    (equivalently the free scalar inverse of the derivative source).\n
    Each order `deriv[mu] = 2 m + e` is split into an even power of
    `d_even(k_mu)` and, for an odd `deriv[mu]` (`e = 1`), one symmetric
    (central) difference factor `i sin(k_mu)`.  `even_deriv_kernel` selects
    `d_even`:\n
    - `"half"` (the default, also `None`): `d_even = -4 sin^2(k_mu / 2)`, minus
      the `mu` term of the `D(k)` used by `free_scalar_mom_invert`, so even
      orders are the laplacian powers;\n
    - `"central"`: `d_even = -sin^2(k_mu)`, the square of `i sin(k_mu)`, so
      every order is a power of the symmetric difference factor.\n
    In particular `deriv=[1, 0, 0, 0]` is the central difference in the `x`
    direction (independent of `even_deriv_kernel`), `deriv=[2, 0, 0, 0]` is
    minus the `x` term of the laplacian with `"half"`, and
    `deriv=[2, 0, 0, 0]` with `"central"` is the square of that central
    difference, `(f(x+2) - 2 f(x) + f(x-2)) / 4`.\n
    `deriv=None` is equivalent to `[0, 0, 0, 0]`, in which case the result
    equals `free_scalar_invert(src, mass, ...)`.
    """
    fft_f = q.mk_fft(is_forward=True, is_normalizing=True, mode_fft=mode_fft)
    fft_b = q.mk_fft(is_forward=False, is_normalizing=True, mode_fft=mode_fft)
    f = fft_f * src
    q.free_scalar_deriv_mom(f, deriv, momtwist, even_deriv_kernel=even_deriv_kernel)
    q.free_scalar_mom_invert(f, mass, momtwist)
    sol = fft_b * f
    return sol
