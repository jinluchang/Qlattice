# cython: binding=True, embedsignature=True

"""
Module ``qlat.vector_utils``
=============================\n
Cython bindings for the legacy vector-utility routines (the former
``qlat/cqlat/vector_utils.cpp`` exports): propagator/gauge-field
comparison, gwu/qlat file I/O, source construction and meson
contraction helpers.\n
Documentation: ``docs/qlat/qlat_vector_utils.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .propagator cimport Prop
from .qcd cimport GaugeField

### -------------------------------------------------------------------

def diff_gauge(GaugeField g0, GaugeField g1):
    return cc.py_diff_gauge(g0.xxx().val(), g1.xxx().val())

def diff_prop(Prop p0, Prop p1):
    cc.py_diff_prop(p0.xxx().val(), p1.xxx().val())

### -------------------------------------------------------------------
### gwu / qlat file I/O

def load_gwu_link(GaugeField g0, path):
    cdef bytes path_b = as_bytes(path)
    cc.py_load_gwu_link(<const char*>path_b, g0.xxx().val(), True)

def save_gwu_prop(Prop prop, path):
    cdef bytes path_b = as_bytes(path)
    cc.py_save_gwu_prop(<const char*>path_b, prop.xxx().val())

def load_gwu_prop(Prop prop, path):
    cdef bytes path_b = as_bytes(path)
    cc.py_load_gwu_prop(<const char*>path_b, prop.xxx().val())

def save_gwu_noiP(Prop prop, path):
    cdef bytes path_b = as_bytes(path)
    cc.py_save_gwu_noiP(<const char*>path_b, prop.xxx().val())

def load_gwu_noiP(Prop prop, path):
    cdef bytes path_b = as_bytes(path)
    cc.py_load_gwu_noiP(<const char*>path_b, prop.xxx().val())

def load_qlat_link(GaugeField g0, path):
    cdef bytes path_b = as_bytes(path)
    cc.py_load_qlat_link(path_b, g0.xxx().val())

def save_qlat_link(GaugeField g0, path):
    cdef bytes path_b = as_bytes(path)
    cc.py_save_qlat_link(path_b, g0.xxx().val())

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
    cdef bytes filename_b = as_bytes(filename)
    cdef bytes info_b = as_bytes(info)
    cc.py_meson_corrE(p0.xxx().val(), p1.xxx().val(), g0, g1, filename_b,
                   mom_c.xx, invmode, tini, info_b, shift_end)

def corr_dat_create(filename, key_T, dimN, info="NONE"):
    cdef bytes filename_b = as_bytes(filename)
    cdef bytes key_T_b = as_bytes(key_T)
    cdef bytes dimN_b = as_bytes(dimN)
    cdef bytes info_b = as_bytes(info)
    cc.corr_dat_create(filename_b, key_T_b, dimN_b, info_b)

def corr_dat_info(filename, info="NONE"):
    cdef bytes filename_b = as_bytes(filename)
    cdef bytes info_b = as_bytes(info)
    cc.corr_dat_info(filename_b, info_b)

### -------------------------------------------------------------------

def prop4d_conj(Prop prop, cc.Int rotate=1):
    cc.py_prop4d_conj(prop.xxx().val(), rotate)

def prop4d_src_gamma(Prop prop, cc.Int g0=0, cc.Int Conj=0):
    cc.py_prop4d_src_gamma(prop.xxx().val(), g0, Conj)

def prop4d_sink_gamma(Prop prop, cc.Int g0=0, cc.Int Conj=0):
    cc.py_prop4d_sink_gamma(prop.xxx().val(), g0, Conj)
