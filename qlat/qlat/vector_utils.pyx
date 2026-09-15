# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.vector_utils``
=============================

Cython bindings for the legacy vector-utility routines (the former
``qlat/cqlat/vector_utils.cpp`` exports): propagator/gauge-field
comparison, gwu/qlat file I/O, source construction and meson
contraction helpers.

Documentation: ``docs/qlat/qlat_vector_utils.md``

.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .propagator cimport Prop
from .qcd cimport GaugeField

import qlat_utils as q

### -------------------------------------------------------------------

def diff_gauge(GaugeField g0, GaugeField g1):
    return cc.py_diff_gauge(g0.xxx().val(), g1.xxx().val())

def diff_prop(Prop p0, Prop p1):
    cc.py_diff_prop(p0.xxx().val(), p1.xxx().val())

### -------------------------------------------------------------------
### gwu / qlat file I/O

def load_gwu_link(GaugeField g0, path):
    cc.py_load_gwu_link(path, g0.xxx().val(), True)

def save_gwu_prop(Prop prop, path):
    cc.py_save_gwu_prop(path, prop.xxx().val())

def load_gwu_prop(Prop prop, path):
    cc.py_load_gwu_prop(path, prop.xxx().val())

def save_gwu_noiP(Prop prop, path):
    cc.py_save_gwu_noiP(path, prop.xxx().val())

def load_gwu_noiP(Prop prop, path):
    cc.py_load_gwu_noiP(path, prop.xxx().val())

def load_qlat_link(GaugeField g0, path):
    cc.py_load_qlat_link(path, g0.xxx().val())

def save_qlat_link(GaugeField g0, path):
    cc.py_save_qlat_link(path, g0.xxx().val())

### -------------------------------------------------------------------
### source construction

def random_point_src(Prop prop, cc.Int seed=0):
    cc.py_random_point_src(prop.xxx().val(), seed)

def make_point_prop(Prop prop, sp=None):
    cdef Coordinate sp_c
    if sp is None:
        sp_c = Coordinate()
    else:
        sp_c = Coordinate(sp)
    cc.py_make_point_prop(prop.xxx().val(), sp_c.xx)

def make_volume_src(Prop prop, cc.Int seed=0, cc.Int mix_color=0,
                    cc.Int mix_spin=0, cc.Int tini=-1):
    cc.py_make_volume_src(prop.xxx().val(), seed, mix_color, mix_spin, tini)

def local_sequential_source(Prop res, Prop src, tseq, cc.Int gammai=-1):
    cdef cc.vector[cc.Int] tseq_v = cc.vector[cc.Int]()
    tseq_v.resize(len(tseq))
    cdef Py_ssize_t i
    for i in range(len(tseq)):
        tseq_v[i] = tseq[i]
    cc.py_local_sequential_source(res.xxx().val(), src.xxx().val(), tseq_v, gammai)

### -------------------------------------------------------------------

def meson_corr(Prop p0, Prop p1, filename, cc.Int g0, cc.Int g1, cc.Int tini=0,
               cc.Int invmode=1, info="NONE", cc.Int shift_end=1, mom=None):
    cdef Coordinate mom_c
    if mom is None:
        mom_c = Coordinate()
    else:
        mom_c = Coordinate(mom)
    cc.py_meson_corrE(p0.xxx().val(), p1.xxx().val(), g0, g1, filename,
                   mom_c.xx, invmode, tini, info, shift_end)

def corr_dat_create(filename, key_T, dimN, info="NONE"):
    cc.corr_dat_create(filename, key_T, dimN, info)

def corr_dat_info(filename, info="NONE"):
    cc.corr_dat_info(filename, info)

### -------------------------------------------------------------------

def prop4d_conj(Prop prop, cc.Int rotate=1):
    cc.py_prop4d_conj(prop.xxx().val(), rotate)

def prop4d_src_gamma(Prop prop, cc.Int g0=0, cc.Int Conj=0):
    cc.py_prop4d_src_gamma(prop.xxx().val(), g0, Conj)

def prop4d_sink_gamma(Prop prop, cc.Int g0=0, cc.Int Conj=0):
    cc.py_prop4d_sink_gamma(prop.xxx().val(), g0, Conj)
