# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

import os

os.environ['CUBACORES'] = '0'

from . cimport everything as cc
from qlat_utils.all cimport *
cimport numpy
import numpy as np
import qlat_utils as q

@q.timer
def test_integration_multi_dimensional():
    cc.test_integrationMultidimensional()

@q.timer
def clear_muon_line_interpolations():
    """
    interface to C++ function `clear_muon_line_interpolations`
    """
    cc.clear_muon_line_interpolations()

@q.timer
def compute_save_muonline_interpolation(cc.std_string path, cc.std_vector[cc.Int] dims, cc.RealD epsabs, cc.RealD epsrel):
    """
    compute and save the muonline interpolation at `path`.
    #
    interface to C++ function `compute_save_muonline_interpolation_cc`
    #
    default: epsabs=1e-8 ; epsrel=1e-3
    dims = [ 6, 6, 6, 6, 6, ] ~ [ 16, 16, 16, 16, 16, ]
    """
    return cc.compute_save_muonline_interpolation_cc(path, dims, epsabs, epsrel)

@q.timer
def load_multiple_muonline_interpolations(cc.std_string path, cc.std_vector[cc.Long] idx_list):
    """
    load muonline interpolation saved at `path`.
    #
    interface to C++ function `load_multiple_muonline_interpolations`
    #
    May load multiple interpolations if under path we have:
    f"{path}/{idx:010d}/checkpoint" where f"{path}/{idx:010d}" is the path for `idx` range from 0 up to some `limit`.
    Load all if idx_list is empty.
    """
    return cc.load_multiple_muonline_interpolations(path, idx_list)

@q.timer
def calc_muon_line_m(CoordinateD x, CoordinateD y, cc.RealD epsabs, cc.RealD epsrel):
    """
    return ret
    ret is 4D np.array with total len=3*4*4*4 and dtype=np.float64
    ret[i, rho, sigma, lambda] = \mathcal M_{i,\\rho,\\sigma,\\lambda}(x,y,z)
    See Eq. (9) of https://arxiv.org/pdf/2304.04423.pdf
    #
    interface to C++ function `muon_line_sym_py`
    default: epsabs=1e-8 ; epsrel=1e-3
    """
    cdef cc.std_vector[cc.RealD] vec = cc.muon_line_sym_py(x.xx, y.xx, epsabs, epsrel)
    cdef numpy.ndarray arr = np.ascontiguousarray(vec, dtype=np.float64).reshape(3, 4, 4, 4)
    return arr

def get_muon_line_m(CoordinateD x, CoordinateD y, CoordinateD z, cc.Int idx, cc.RealD epsabs, cc.RealD epsrel):
    """
    return ret
    ret is 4D np.array with total len=3*4*4*4 and dtype=np.float64
    ret[i, rho, sigma, lambda] = \mathcal M_{i,\\rho,\\sigma,\\lambda}(x,y,z)
    See Eq. (9) of https://arxiv.org/pdf/2304.04423.pdf
    #
    interface to C++ function `get_muon_line_m_py`
    if idx < 0: calculate instead of using loaded interpolation.
    else using the muonline interpolation loaded as `idx` (epsabs and epsrel are ignored).
    """
    cdef cc.std_vector[cc.RealD] vec = cc.get_muon_line_m_py(x.xx, y.xx, z.xx, idx, epsabs, epsrel)
    cdef numpy.ndarray arr = np.ascontiguousarray(vec, dtype=np.float64).reshape(3, 4, 4, 4)
    return arr

def get_muon_line_m_extra(CoordinateD x, CoordinateD y, CoordinateD z, cc.Int tag):
    """
    return ret
    ret is 4D np.array with total len=3*4*4*4 and dtype=np.float64
    ret[i, rho, sigma, lambda] = \mathcal M_{i,\\rho,\\sigma,\\lambda}(x,y,z)
    See Eq. (9) of https://arxiv.org/pdf/2304.04423.pdf
    #
    interface to C++ function `get_muon_line_m_extra_py`
    tag = 0 sub
    tag = 1 nosub
    """
    cdef cc.std_vector[cc.RealD] vec = cc.get_muon_line_m_extra_py(x.xx, y.xx, z.xx, tag)
    cdef numpy.ndarray arr = np.ascontiguousarray(vec, dtype=np.float64).reshape(3, 4, 4, 4)
    return arr

def get_muon_line_m_extra_lat(Coordinate x, Coordinate y, Coordinate z, Coordinate total_site, cc.RealD a, cc.Int tag):
    """
    return ret
    ret is 4D np.array with total len=3*4*4*4 and dtype=np.float64
    ret[i, rho, sigma, lambda] = \mathcal M_{i,\\rho,\\sigma,\\lambda}(x,y,z)
    See Eq. (9) of https://arxiv.org/pdf/2304.04423.pdf
    #
    interface to C++ function `get_muon_line_m_extra_lat_py`
    a is muon_mass in lattice unit.
    tag = 0 sub
    tag = 1 nosub
    """
    cdef cc.std_vector[cc.RealD] vec = cc.get_muon_line_m_extra_lat_py(x.xx, y.xx, z.xx, total_site.xx, a, tag)
    cdef numpy.ndarray arr = np.ascontiguousarray(vec, dtype=np.float64).reshape(3, 4, 4, 4)
    return arr
