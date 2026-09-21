# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.field_utils``
============================\n
Field utility functions: communication plans, coordinate shifting, and
field shuffling.  The pure-Python helpers (field expansion, halo refresh,
FFT, norms, and element-wise square root) live in
``qlat.field_utils_utils``.\n
Documentation: ``docs/qlat/qlat_field_utils.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .field_base cimport FieldBase
from .field_types cimport (
        FieldInt8t,
        FieldRealD,
        FieldComplexD,
        FieldChar,
        )
from .selected_field_types cimport SelectedFieldRealD
from .selected_points_types cimport SelectedPointsRealD

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

class q:
    from qlat_utils import (
        timer,
    )
    from .field_base_utils import (
        Field,
    )

cdef class CommMarks(FieldInt8t):

    def __init__(self, Geometry geo=None, int multiplicity=0):
        super().__init__(geo, multiplicity)

###

cdef class CommPlan:

    def __imatmul__(self, CommPlan v1):
        self.xx = v1.xx
        return self

    def copy(self, is_copying_data=True):
        x = type(self)()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

###

def make_field_expand_comm_plan(CommMarks comm_marks):
    """
    comm_marks is of type Field(ElemTypeInt8t)
    """
    cdef CommPlan cp = CommPlan()
    cp.xx = cc.py_make_comm_plan(comm_marks.xx)
    return cp

def set_marks_field_all(CommMarks comm_marks, Geometry geo, int multiplicity,
                        tag):
    """
    Set the standard expansion marks for ``comm_marks``.
    """
    cc.py_set_marks_field_all(comm_marks.xx, geo.xx, multiplicity, tag)

def mk_phase_field(Geometry geo, lmom):
    """
    lmom is in lattice momentum unit
    exp(i * 2*pi/L * lmom \cdot xg )
    """
    cdef CoordinateD lmom_d = CoordinateD(lmom)
    cdef FieldComplexD f = q.Field(ElemTypeComplexD, geo, 1)
    cc.py_set_phase_field(f.xx, lmom_d.xx)
    return f

@q.timer
def sqrt_selected_points_real_d(SelectedPointsRealD f):
    cdef SelectedPointsRealD f_ret = f.copy(is_copying_data=False)
    cc.set_sqrt_field(f_ret.xx, f.xx)
    return f_ret

@q.timer
def sqrt_selected_field_real_d(SelectedFieldRealD f):
    cdef SelectedFieldRealD f_ret = f.copy(is_copying_data=False)
    cc.set_sqrt_field(f_ret.xx, f.xx)
    return f_ret

@q.timer
def sqrt_field_real_d(FieldRealD f):
    cdef FieldRealD f_ret = f.copy(is_copying_data=False)
    cc.set_sqrt_field(f_ret.xx, f.xx)
    return f_ret

@q.timer
def field_char_shift(FieldChar f, Coordinate shift):
    """
    return `sf`
    `sf` is the new shifted Field.
    `f` is not changed.
    `shift` is the coordinate to shift the field (or None, which means no shift)
    roughly shifted_field[(xg + shift) % total_site] == original_field[xg]
    """
    cdef FieldChar sf = FieldChar()
    cc.field_shift(sf.xx, f.xx, shift.xx)
    return sf

@q.timer
def field_shift(FieldBase f, Coordinate shift):
    cdef FieldChar fc = FieldChar()
    f.swap_cast(fc)
    cdef FieldChar fc_shifted = field_char_shift(fc, shift)
    sf = type(f)()
    sf.swap_cast(fc_shifted)
    f.swap_cast(fc)
    return sf

@q.timer
def shuffle_field_char(FieldChar f, Coordinate new_size_node):
    """
    return f_list
    """
    cdef cc.std_vector[cc.Field[cc.Char]] fs
    fs = cc.std_vector[cc.Field[cc.Char]]()
    cc.shuffle_field(fs, f.xx, new_size_node.xx)
    cdef cc.Int num_field = <cc.Int>fs.size()
    cdef cc.Int i
    cdef FieldChar f1
    cdef list f_list = []
    for i in range(num_field):
        f1 = FieldChar()
        cc.qswap(f1.xx, fs[i])
        f_list.append(f1)
    return f_list

@q.timer
def shuffle_field_char_back(FieldChar fc, list f_list, Coordinate new_size_node):
    """
    Modify `fc` in place.
    NOTE: `fc` needs to have correct size.
    """
    cdef cc.std_vector[cc.Field[cc.Char]] fs
    cdef cc.Int num_field = <cc.Int>len(f_list)
    cdef cc.Int i
    cdef FieldChar f1
    fs = cc.std_vector[cc.Field[cc.Char]](num_field)
    for i in range(num_field):
        f1 = f_list[i]
        cc.qswap(f1.xx, fs[i])
    cc.shuffle_field_back(fc.xx, fs, new_size_node.xx)
    for i in range(num_field):
        f1 = f_list[i]
        cc.qswap(f1.xx, fs[i])

@q.timer
def shuffle_field(FieldBase f, Coordinate new_size_node):
    cdef FieldChar fc = FieldChar()
    f.swap_cast(fc)
    cdef list fc_list = shuffle_field_char(fc, new_size_node)
    f.swap_cast(fc)
    cdef list f_list = []
    for i in range(len(fc_list)):
        f1 = type(f)()
        f1.swap_cast(fc_list[i])
        f_list.append(f1)
    return f_list

@q.timer
def shuffle_field_back(FieldBase f, list f_list, Coordinate new_size_node):
    """
    Modify `f` in place
    NOTE: `f` needs to have correct size.
    """
    cdef FieldChar f1
    cdef list fc_list = []
    for i in range(len(f_list)):
        f1 = FieldChar()
        f_list[i].swap_cast(f1)
        fc_list.append(f1)
    cdef FieldChar fc = FieldChar()
    f.swap_cast(fc)
    shuffle_field_char_back(fc, fc_list, new_size_node)
    f.swap_cast(fc)
    for i in range(len(fc_list)):
        f_list[i].swap_cast(fc_list[i])