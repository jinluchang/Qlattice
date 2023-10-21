# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .qcd cimport GaugeField
from .propagator cimport Prop

import cqlat as c
import qlat_utils as q
import numpy as np

def gf_ape_smear(GaugeField gf, double alpha, int steps=1):
    return cc.gf_ape_smear(gf.xxx().val(), gf.xxx().val(), alpha, steps)

def gf_spatial_ape_smear(GaugeField gf, double alpha, int steps=1):
    """
    used value: alpha = 0.5, steps = 30
    """
    return c.gf_spatial_ape_smear(gf, gf, alpha, steps)

def gf_hyp_smear(GaugeField gf, double alpha1, double alpha2, double alpha3):
    """
    values in paper is 0.75 0.6 0.3
    10.1103/PhysRevD.64.034504 Eq(4)
    """
    return c.gf_hyp_smear(gf, gf, alpha1, alpha2, alpha3)

def prop_smear(Prop prop, GaugeField gf1,
               double coef, int step, mom=None,
               cc.bool smear_in_time_dir=False,
               int mode_smear=1):
    """
    gf_ape = gf_spatial_ape_smear(gf, 0.5, 30)
    gf1 = mk_left_expanded_gauge_field(gf_ape)
    mom: momentum smearing in lattice unit 1/a (not unit of lattice momentum 2 pi / L / a)
    smear params:
    24D: coef = 0.9375, step = 10
    48I: coef = 0.9375, step = 29
    mom = 0.5 * mom of the corresponding hadron
    """
    if mom is None:
        mom = [ 0.0, 0.0, 0.0, 0.0, ]
    return c.prop_smear(prop, gf1, coef, step, mom, smear_in_time_dir, mode_smear)
