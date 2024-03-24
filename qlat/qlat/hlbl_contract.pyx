# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
from qlat_utils.all cimport *
cimport numpy
import numpy as np
import qlat_utils as q
from .field_selection cimport FieldSelection, PointsSelection
from .field_types cimport FieldRealD, FieldComplexD
from .selected_field_types cimport SelectedFieldRealD, SelectedFieldComplexD
from .selected_points_types cimport SelectedPointsRealD, SelectedPointsComplexD
from .propagator cimport SelProp
from .geometry cimport Geometry

cdef numpy.ndarray sl_arr_from_sl_table(const cc.SlTable& x):
    """
    return sl_arr
    sl_arr.shape == (s_limit, l_limit,)
    sl_arr.dtype == np.complex128
    """
    cdef cc.Long s_limit = x.s_limit
    cdef cc.Long l_limit = x.l_limit
    cdef numpy.ndarray sl_arr = np.zeros((s_limit, l_limit,), dtype=np.complex128)
    cdef cc.PyComplexD[:] ptr = sl_arr.ravel()
    cdef cc.Long table_size = x.table.size()
    cc.memcpy(&ptr[0], x.table.data(), table_size * sizeof(cc.PyComplexD))
    return sl_arr

@q.timer
def mk_m_z_field_tag(
        FieldSelection fsel,
        Coordinate xg_x,
        Coordinate xg_y,
        const cc.RealD a,
        const cc.Int tag
        ):
    """
    return smf_d
    smf is SelectedField("double") as the muon-line-field of vertex z
    a is the lattice spacing (when muon_mass = 1). In lattice unit one should use muon_mass * a for this parameter
    tag = 0: sub
    tag = 1: nosub
    """
    cdef SelectedFieldRealD smf_d = SelectedFieldRealD(fsel)
    cc.set_m_z_field_tag(smf_d.xx, fsel.xx, xg_x.xx, xg_y.xx, a, tag)
    return smf_d

def contract_four_pair_labels(list tags):
    return cc.contract_four_pair_labels(tags)

def contract_two_plus_two_pair_labels():
    return cc.contract_two_plus_two_pair_labels()

@q.timer
def contract_four_pair_no_glb_sum(
        cc.PyComplexD coef,
        SelectedPointsRealD psel_prob,
        SelectedFieldRealD fsel_prob,
        const cc.Long idx_xg_x,
        const cc.Long idx_xg_y,
        SelectedFieldRealD smf_d,
        SelProp sprop_x,
        SelProp sprop_y,
        const cc.Int inv_type,
        const cc.std_vector[cc.std_string]& tags,
        const cc.Long r_sq_limit,
        const cc.RealD muon_mass,
        const cc.RealD z_v
        ):
    """
    return lsl_arr
    #
    lsl_arr.shape == (len(labels), s_limit, l_limit,)
    labels = contract_four_pair_labels(tags)
    #
    default coef = 1.0
    inv_type = 0 : light quark
    inv_type = 1 : strange quark
    tags can be [ "ref-far", "ref-center", "ref-close", ]
    #
    glb_sum for SlTable not yet performed
    #
    x, y are two sampled points (elem of `psel`). `psel_prob` factors are already included.
    z, x_op are summed over with in `fsel`. `fsel_prob` factors are already included.
    """
    cdef PointsSelection psel = psel_prob.psel
    cdef FieldSelection fsel = fsel_prob.fsel
    assert fsel is smf_d.fsel
    assert fsel is sprop_x.fsel
    assert fsel is sprop_y.fsel
    cdef cc.std_vector[cc.SlTable] sl_table_vec = cc.contract_four_pair_no_glb_sum(
            cc.ccpy_d(coef),
            psel.xx,
            psel_prob.xx,
            fsel.xx,
            fsel_prob.xx,
            idx_xg_x,
            idx_xg_y,
            smf_d.xx,
            sprop_x.xx,
            sprop_y.xx,
            inv_type,
            tags,
            r_sq_limit,
            muon_mass,
            z_v,
            )
    cdef cc.Long s = sl_table_vec.size()
    cdef list sl_arr_list = []
    for i in range(s):
        sl_arr = sl_arr_from_sl_table(sl_table_vec[i])
        sl_arr_list.append(sl_arr)
    return np.array(sl_arr_list, dtype=np.complex128)

@q.timer
def contract_two_plus_two_pair_no_glb_sum(
        cc.PyComplexD coef,
        SelectedPointsRealD psel_prob,
        SelectedPointsRealD psel_lps_prob,
        const cc.Long idx_xg_x,
        SelectedPointsComplexD lps_hvp_x,
        SelectedPointsComplexD edl_list_c,
        const cc.Long r_sq_limit,
        const cc.RealD muon_mass,
        const cc.RealD z_v,
        ):
    """
    return n_points_in_r_sq_limit, n_points_computed, lsl_arr
    #
    lsl_arr.shape == (len(labels), s_limit, l_limit,)
    labels = contract_two_plus_two_pair_labels()
    #
    hvp point source at x (rho)
    hvp point sink at y (sigma)
    hvp with external loop source at z (lambda)
    hvp with external loop sink summed over (i)
    psel_edl[k] = xg_z
    edl_list[k][i * 4 + lambda]
    -hvp_x.get_elem(xl_y, sigma * 4 + rho)
    #
    glb_sum for SlTable not yet performed
    #
    x is sampled point (elem of `psel`). `psel_prob` factor is already included.
    y is summed over with adaptive sampling. Weight factors is already included.
    z is sampled point (elem of `psel`) and summed over. `psel_prob` factor is already included.
    x_op is summed over all points.
    """
    cdef Geometry geo = psel_prob.psel.geo
    cdef PointsSelection psel = psel_prob.psel
    cdef PointsSelection psel_lps = psel_lps_prob.psel
    cdef cc.Long n_points_in_r_sq_limit = 0
    cdef cc.Long n_points_computed = 0
    cdef cc.std_vector[cc.SlTable] sl_table_vec = cc.contract_two_plus_two_pair_no_glb_sum(
            n_points_in_r_sq_limit,
            n_points_computed,
            cc.ccpy_d(coef),
            geo.xx,
            psel.xx,
            psel_prob.xx,
            psel_lps.xx,
            psel_lps_prob.xx,
            idx_xg_x,
            lps_hvp_x.xx,
            edl_list_c.xx,
            r_sq_limit,
            muon_mass,
            z_v,
            )
    cdef cc.Long s = sl_table_vec.size()
    cdef list sl_arr_list = []
    for i in range(s):
        sl_arr = sl_arr_from_sl_table(sl_table_vec[i])
        sl_arr_list.append(sl_arr)
    lsl_arr = np.array(sl_arr_list, dtype=np.complex128)
    return n_points_in_r_sq_limit, n_points_computed, lsl_arr
