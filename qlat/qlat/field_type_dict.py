"""
Module ``qlat.field_type_dict``
================================\n
Registry of field element-type classifications used throughout Qlattice.\n
Defines lookup dictionaries (``field_type_dict``, ``selected_field_type_dict``,
``selected_points_type_dict``) and pre-built lists of element-type categories
(``field_ctypes_complex``, ``field_ctypes_double``, ``field_ctypes_float``,
``field_ctypes_long``, ``field_ctypes_char``) that other modules use to
determine which C++ code path to invoke for a given field element type.\n
Documentation: ``docs/qlat/qlat_field_type_dict.md``\n
.. note:: Update the documentation when updating this source file.
"""

__all__ = [
    "field_type_dict",
    "selected_field_type_dict",
    "selected_points_type_dict",
    "field_ctypes_complex",
    "field_ctypes_double",
    "field_ctypes_float",
    "field_ctypes_long",
    "field_ctypes_char",
]

### -------------------------------------------------------------------

class q:
    from qlat_utils import (
        ElemTypeChar,
        ElemTypeColorMatrix,
        ElemTypeComplexD,
        ElemTypeComplexF,
        ElemTypeInt64t,
        ElemTypeInt8t,
        ElemTypeIsospinMatrix,
        ElemTypeLong,
        ElemTypeNonRelWilsonMatrix,
        ElemTypeRealD,
        ElemTypeRealF,
        ElemTypeSpinMatrix,
        ElemTypeWilsonMatrix,
        ElemTypeWilsonVector,
    )

### -------------------------------------------------------------------

field_type_dict = {}

selected_field_type_dict = {}

selected_points_type_dict = {}

### -------------------------------------------------------------------

field_ctypes_complex = [
    q.ElemTypeColorMatrix,
    q.ElemTypeWilsonMatrix,
    q.ElemTypeNonRelWilsonMatrix,
    q.ElemTypeIsospinMatrix,
    q.ElemTypeSpinMatrix,
    q.ElemTypeWilsonVector,
    q.ElemTypeComplexD,
]

field_ctypes_complex_f = [
    q.ElemTypeComplexF,
]

field_ctypes_complex_or_complex_f = field_ctypes_complex + field_ctypes_complex_f

field_ctypes_double = [
    q.ElemTypeColorMatrix,
    q.ElemTypeWilsonMatrix,
    q.ElemTypeNonRelWilsonMatrix,
    q.ElemTypeIsospinMatrix,
    q.ElemTypeSpinMatrix,
    q.ElemTypeWilsonVector,
    q.ElemTypeComplexD,
    q.ElemTypeRealD,
]

field_ctypes_float = [
    q.ElemTypeComplexF,
    q.ElemTypeRealF,
]

field_ctypes_double_or_float = field_ctypes_double + field_ctypes_float

field_ctypes_long = [
    q.ElemTypeLong,
    q.ElemTypeInt64t,
]

field_ctypes_char = [
    q.ElemTypeChar,
    q.ElemTypeInt8t,
]
