# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.field_double``
==============================\n
Element-wise operations on real-double (``ElemTypeRealD``) lattice
fields, including type conversion, comparison, inversion, and
multiplication.\n
Documentation: ``docs/qlat/qlat_field_double.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_base cimport FieldBase
from .field_types cimport FieldRealD
from .field_types cimport FieldComplexD

import cqlat as c

def set_double_from_complex(field, FieldComplexD cf):
    assert isinstance(field, FieldBase)
    assert field.ctype == ElemTypeRealD
    assert cf.ctype == ElemTypeComplexD
    field._cc_set_double_from_complex(cf)

def set_complex_from_double(field, FieldRealD sf):
    assert isinstance(field, FieldBase)
    assert field.ctype == ElemTypeComplexD
    assert sf.ctype == ElemTypeRealD
    field._cc_set_complex_from_double(sf)

def set_abs_from_complex(field, FieldComplexD cf):
    assert isinstance(field, FieldBase)
    assert field.ctype == ElemTypeRealD
    assert cf.ctype == ElemTypeComplexD
    field._cc_set_abs_from_complex(cf)

def set_ratio_double(field, FieldRealD sf1, FieldRealD sf2):
    assert isinstance(field, FieldBase)
    assert field.ctype == ElemTypeRealD
    field._cc_set_ratio_double(sf1, sf2)

def less_than_double(field, FieldRealD sf2, FieldRealD mask):
    assert isinstance(field, FieldBase)
    assert field.ctype == ElemTypeRealD
    field._cc_less_than_double(sf2, mask)

def invert_double(field):
    assert field.ctype == ElemTypeRealD
    c.invert_double_field(field)

def multiply_double(field, factor):
    assert isinstance(field, FieldBase)
    assert isinstance(factor, FieldBase)
    assert factor.ctype is ElemTypeRealD
    c.multiply_double_field(field, factor)
