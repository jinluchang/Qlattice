# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.smear``
======================\n
Gauge-field smearing (APE, HYP) and propagator/fermion-field smearing
operations for lattice QCD calculations.\n
Documentation: ``docs/qlat/qlat_smear.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .qcd cimport GaugeField
from .propagator cimport Prop, FermionField4d

class q:
    from qlat_utils import (
        timer,
        CoordinateD,
    )

import numpy as np

def gf_ape_smear(GaugeField gf, double alpha, cc.Long num_step=1):
    cdef GaugeField gf1 = gf.copy()
    cc.gf_ape_smear(gf1.xxx().val(), gf.xxx().val(), alpha, num_step)
    return gf1

def gf_ape_smear_no_proj(GaugeField gf, cc.RealD alpha, cc.RealD accumulated_alpha):
    """
    Similar to `gf_ape_smear` but with the SU(3) projection step omitted.
    `accumulated_alpha` is the target total effective smearing; every step
    updates the links as `U' = U + (a / 6) * S` (no `(1 - a)` factor) with
    `a = alpha * 3 / M`, where `M` is the global maximum of
    `Tr(U_mu(x)^dag U_mu(x))` over all links, so that `alpha` is the
    effective smearing step of the largest link. The field grows under this
    kernel (`U' ~ (1 + a) U` for near-unitary links), so `M > 3` and hence
    `a <= alpha`. The weight decays geometrically as the amplitude grows, so
    the accumulated smearing converges; a target above that limit is never
    reached and the loop stops on the no-progress guard, leaving a very large,
    non-unitary field. The field is never rescaled; the last step is truncated
    so that `sum(a)` does not exceed `accumulated_alpha`.
    `alpha = 0` or `accumulated_alpha = 0` returns the input unchanged.
    Returns a new gauge field; `gf` is not modified.
    """
    cdef GaugeField gf1 = gf.copy()
    cc.gf_ape_smear_no_proj(gf1.xxx().val(), gf.xxx().val(), alpha, accumulated_alpha)
    return gf1

def gf_ape_smear_no_proj_rescaling(GaugeField gf, double alpha, cc.Long num_step=1):
    """
    Similar to `gf_ape_smear` but with the SU(3) projection step omitted.
    Each step uses the convex kernel `U' = (1 - alpha) U + (alpha / 6) * S`,
    then multiplies the whole field by the real scalar
    `s = sqrt(3 / <N^2>)`, where `<N^2>` is the global (MPI-reduced) average
    of `Tr(U_mu(x)^dag U_mu(x))` over all links, so that `<N^2> = 3` (the
    value for a unitary link) after every step, including the last.
    Because the whole field is rescaled, `alpha = 0` only removes the staple
    mixing; the rescale still applies, so it is not the identity for a
    non-unitary field. `num_step = 0` returns an unchanged copy.
    Returns a new gauge field; `gf` is not modified.
    """
    cdef GaugeField gf1 = gf.copy()
    cc.gf_ape_smear_no_proj_rescaling(gf1.xxx().val(), gf.xxx().val(), alpha, num_step)
    return gf1

def gf_ape_smear_spatial(GaugeField gf, double alpha, cc.Long num_step=1):
    """
    used value: alpha = 0.5, num_step = 30
    """
    cdef GaugeField gf1 = gf.copy()
    cc.gf_ape_smear_spatial(gf1.xxx().val(), gf.xxx().val(), alpha, num_step)
    return gf1

def gf_hyp_smear(GaugeField gf, double alpha1, double alpha2, double alpha3, cc.Long num_step=1):
    """
    values in paper is 0.75 0.6 0.3
    10.1103/PhysRevD.64.034504 Eq(4)
    """
    cdef GaugeField gf1 = gf.copy()
    cc.gf_hyp_smear(gf1.xxx().val(), gf.xxx().val(), alpha1, alpha2, alpha3, num_step)
    return gf1

@q.timer
def prop_smear(Prop prop, GaugeField gf,
               double coef, cc.Int num_step=1, CoordinateD mom=None,
               cc.Bool smear_in_time_dir=False,
               str mode="qlat"):
    """
    gf_ape = gf_ape_smear_spatial(gf, 0.5, 30)
    gf = mk_left_expanded_field(gf_ape)
    mom: momentum smearing in lattice unit 1/a (not unit of lattice momentum 2 pi / L / a)
    mode: "qlat" (default, qlat convention) or "standard" (legacy `prop_smear`)
    smear params:
    24D: coef = 0.9375, num_step = 10
    48I: coef = 0.9375, num_step = 29
    mom = 0.5 * mom of the corresponding hadron
    """
    cdef Prop prop1 = prop.copy()
    if mom is None:
        mom = CoordinateD()
    if mode == "standard":
        cc.prop_smear(prop1.xxx().val(), gf.xxx().val(),
                      coef, num_step, mom.xx, smear_in_time_dir)
    elif mode == "qlat":
        cc.prop_smear_qlat_convention(prop1.xxx().val(), gf.xxx().val(),
                                      coef, num_step, mom.xx, smear_in_time_dir, 1)
    else:
        raise Exception(f"prop_smear(prop, gf, coef={coef}, num_step={num_step}, mom={mom}, smear_in_time_dir={smear_in_time_dir}, mode={mode})")
    return prop1

@q.timer
def prop_spatial_smear_no_comm(list ff_list, GaugeField gf, cc.RealD coef, cc.Long step, CoordinateD mom=None):
    """
    Internal low-level kernel.
    `gf` and each of `ff_list` should contain entire time slices.
    No communication will be performed.
    """
    if mom is None:
        mom = q.CoordinateD()
    cdef cc.Int num_field = len(ff_list)
    cdef cc.std_vector[cc.FermionField4d] ff_vec
    cdef FermionField4d ff
    ff_vec.resize(num_field)
    for i in range(num_field):
        assert isinstance(ff_list[i], FermionField4d)
        ff = ff_list[i]
        cc.qswap(ff.xx, ff_vec[i])
    cc.prop_spatial_smear_no_comm(ff_vec, gf.xxx().val(), coef, step, mom.xx)
    for i in range(num_field):
        ff = ff_list[i]
        cc.qswap(ff.xx, ff_vec[i])

