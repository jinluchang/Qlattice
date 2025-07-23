__all__ = [
        'field_type_dict',
        'selected_field_type_dict',
        'selected_points_type_dict',
        'field_ctypes_complex',
        'field_ctypes_double',
        'field_ctypes_float',
        'field_ctypes_long',
        'field_ctypes_char',
        ]

### -------------------------------------------------------------------

from qlat_utils import *

### -------------------------------------------------------------------

field_type_dict = {}

selected_field_type_dict = {}

selected_points_type_dict = {}

### -------------------------------------------------------------------

field_ctypes_complex = [
        ElemTypeColorMatrix,
        ElemTypeWilsonMatrix,
        ElemTypeNonRelWilsonMatrix,
        ElemTypeIsospinMatrix,
        ElemTypeSpinMatrix,
        ElemTypeWilsonVector,
        ElemTypeComplexD,
        ]

field_ctypes_complex_f = [
        ElemTypeComplexF,
        ]

field_ctypes_complex_or_complex_f = field_ctypes_complex + field_ctypes_complex_f

field_ctypes_double = [
        ElemTypeColorMatrix,
        ElemTypeWilsonMatrix,
        ElemTypeNonRelWilsonMatrix,
        ElemTypeIsospinMatrix,
        ElemTypeSpinMatrix,
        ElemTypeWilsonVector,
        ElemTypeComplexD,
        ElemTypeRealD,
        ]

field_ctypes_float = [
        ElemTypeComplexF,
        ElemTypeRealF,
        ]

field_ctypes_double_or_float = field_ctypes_double + field_ctypes_float

field_ctypes_long = [
        ElemTypeLong,
        ElemTypeInt64t,
        ]

field_ctypes_char = [
        ElemTypeChar,
        ElemTypeInt8t,
        ]
