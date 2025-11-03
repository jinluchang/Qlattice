# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .qcd cimport GaugeField
from .propagator cimport Prop, FermionField4d

import cqlat as c
import qlat_utils as q
import numpy as np

def gf_ape_smear(GaugeField gf, double alpha, int steps=1):
    cdef GaugeField gf1 = gf.copy()
    cc.gf_ape_smear(gf1.xxx().val(), gf.xxx().val(), alpha, steps)
    return gf1

def gf_spatial_ape_smear(GaugeField gf, double alpha, int steps=1):
    """
    used value: alpha = 0.5, steps = 30
    """
    cdef GaugeField gf1 = gf.copy()
    cc.gf_spatial_ape_smear(gf1.xxx().val(), gf.xxx().val(), alpha, steps)
    return gf1

def gf_hyp_smear(GaugeField gf, double alpha1, double alpha2, double alpha3):
    """
    values in paper is 0.75 0.6 0.3
    10.1103/PhysRevD.64.034504 Eq(4)
    """
    cdef GaugeField gf1 = gf.copy()
    cc.gf_hyp_smear(gf1.xxx().val(), gf.xxx().val(), alpha1, alpha2, alpha3)
    return gf1

@q.timer
def prop_smear(Prop prop, GaugeField gf1,
               double coef, int step, CoordinateD mom=None,
               cc.Bool smear_in_time_dir=False,
               int mode_smear=1):
    """
    gf_ape = gf_spatial_ape_smear(gf, 0.5, 30)
    gf1 = mk_left_expanded_field(gf_ape)
    mom: momentum smearing in lattice unit 1/a (not unit of lattice momentum 2 pi / L / a)
    smear params:
    24D: coef = 0.9375, step = 10
    48I: coef = 0.9375, step = 29
    mom = 0.5 * mom of the corresponding hadron
    """
    cdef Prop prop1 = prop.copy()
    if mom is None:
        mom = CoordinateD()
    if mode_smear == 0:
        cc.prop_smear(prop1.xxx().val(), gf1.xxx().val(),
                      coef, step, mom.xx, smear_in_time_dir)
    elif mode_smear >= 1:
        cc.prop_smear_qlat_convension(prop1.xxx().val(), gf1.xxx().val(),
                                      coef, step, mom.xx, smear_in_time_dir, mode_smear)
    else:
        raise Exception(f"prop_smear(prop1, gf1, coef={coef}, step={step}, mom={mom}, smear_in_time_dir={smear_in_time_dir}, mode_smear={mode_smear})")
    return prop1

@q.timer
def prop_spatial_smear_no_comm(list ff_list, GaugeField gf, cc.RealD coef, cc.Long step, CoordinateD mom=None):
    """
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

