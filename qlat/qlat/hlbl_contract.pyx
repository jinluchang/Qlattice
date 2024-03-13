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

cdef numpy.ndarray sl_arr_from_sl_table(const cc.SlTable& x):
    """
    return sl_arr
    sl_arr.shape == (s_limit, l_limit,)
    """
    cdef cc.Long s_limit = x.s_limit
    cdef cc.Long l_limit = x.l_limit
    cdef numpy.ndarray sl_arr = np.zeros((s_limit, l_limit,), dtype=np.complex128)
    cdef cc.PyComplexD* ptr = <cc.PyComplexD*>&x.table[0]
    cdef const cc.PyComplexD[::1] table_view = <cc.PyComplexD[:x.table.size()]>ptr
    sl_arr.ravel()[:] = table_view[:]
    return sl_arr

def set_m_z_field_tag(
        FieldSelection fsel,
        Coordinate xg_x,
        Coordinate xg_y,
        const cc.RealD a,
        const cc.Int tag
        ):
    cdef SelectedFieldRealD smf_d = SelectedFieldRealD(fsel)
    cc.set_m_z_field_tag(smf_d.xx, fsel.xx, xg_x.xx, xg_y.xx, a, tag)
    return smf_d

def contract_four_pair_labels(list tags):
    return cc.contract_four_pair_labels(tags)

def contract_two_plus_two_pair_labels():
    return cc.contract_two_plus_two_pair_labels()

@q.timer
def contract_four_pair(
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
    return [ (label, sl_arr,), ... ]
    #
    default coef = 1.0
    inv_type = 0 : light quark
    inv_type = 1 : strange quark
    tags can be [ "ref-far", "ref-center", "ref-close", ]
    """
    cdef PointsSelection psel = psel_prob.psel
    cdef FieldSelection fsel = fsel_prob.fsel
    assert fsel is smf_d.fsel
    assert fsel is sprop_x.fsel
    assert fsel is sprop_y.fsel
    cdef cc.std_vector[cc.SlTable] sl_table_vec = cc.contract_four_pair(
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
    cdef list labels = contract_four_pair_labels(tags)
    cdef cc.Long s = sl_table_vec.size()
    assert s == len(labels)
    ret = []
    for i, label in enumerate(labels):
        sl_arr = sl_arr_from_sl_table(sl_table_vec[i])
        ret.append((label, sl_arr,))
    return ret

@q.timer
def contract_two_plus_two_pair_no_glb_sum(
        cc.PyComplexD coef,
        SelectedPointsRealD psel_prob,
        FieldRealD rand_prob_sel_field,
        const cc.RealD hvp_sel_threshold,
        const cc.Long idx_xg_x,
        FieldComplexD hvp_x,
        SelectedPointsComplexD edl_list_c,
        const cc.Long r_sq_limit,
        const cc.RealD muon_mass,
        const cc.RealD z_v,
        ):
    """
    hvp point source at x (rho)
    hvp point sink at y (sigma)
    hvp with external loop source at z (lambda)
    hvp with external loop sink summed over (i)
    psel_edl[k] = xg_z
    edl_list[k][i * 4 + lambda]
    -hvp_x.get_elem(xl_y, sigma * 4 + rho)
    #
    glb_sum for SlTable not yet performed
    """
    cdef PointsSelection psel = psel_prob.psel
    cdef cc.Long n_points_in_r_sq_limit = 0
    cdef cc.Long n_points_computed = 0
    cdef cc.std_vector[cc.SlTable] sl_table_vec = cc.contract_two_plus_two_pair_no_glb_sum(
            n_points_in_r_sq_limit,
            n_points_computed,
            cc.ccpy_d(coef),
            psel.xx,
            psel_prob.xx,
            rand_prob_sel_field.xx,
            hvp_sel_threshold,
            idx_xg_x,
            hvp_x.xx,
            edl_list_c.xx,
            r_sq_limit,
            muon_mass,
            z_v,
            )
    cdef list labels = contract_two_plus_two_pair_labels()
    cdef cc.Long s = sl_table_vec.size()
    assert s == len(labels)
    ret = []
    for i, label in enumerate(labels):
        sl_arr = sl_arr_from_sl_table(sl_table_vec[i])
        ret.append((label, sl_arr,))
    return ret
