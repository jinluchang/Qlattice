"""
Module ``qlat.field_double``
==============================\n
Element-wise operations on real-double (``ElemTypeRealD``) lattice
fields, including type conversion, comparison, inversion, and
multiplication.\n
Documentation: ``docs/qlat/qlat_field_double.md``\n
.. note:: Update the documentation when updating this source file.
"""

class q:
    from qlat_utils import (
        ElemTypeComplexD,
        ElemTypeRealD,
    )
    from .field_base import (
        FieldBase,
    )
    from .field_types import (
        FieldComplexD,
        FieldRealD,
    )

def set_double_from_complex(field, cf):
    assert isinstance(field, q.FieldBase)
    assert field.ctype == q.ElemTypeRealD
    assert isinstance(cf, q.FieldComplexD)
    assert cf.ctype == q.ElemTypeComplexD
    field._cc_set_double_from_complex(cf)

def set_complex_from_double(field, sf):
    assert isinstance(field, q.FieldBase)
    assert field.ctype == q.ElemTypeComplexD
    assert isinstance(sf, q.FieldRealD)
    assert sf.ctype == q.ElemTypeRealD
    field._cc_set_complex_from_double(sf)

def set_abs_from_complex(field, cf):
    assert isinstance(field, q.FieldBase)
    assert field.ctype == q.ElemTypeRealD
    assert isinstance(cf, q.FieldComplexD)
    assert cf.ctype == q.ElemTypeComplexD
    field._cc_set_abs_from_complex(cf)

def set_ratio_double(field, sf1, sf2):
    assert isinstance(field, q.FieldBase)
    assert field.ctype == q.ElemTypeRealD
    assert isinstance(sf1, q.FieldRealD)
    assert isinstance(sf2, q.FieldRealD)
    field._cc_set_ratio_double(sf1, sf2)

def less_than_double(field, sf2, mask):
    assert isinstance(field, q.FieldBase)
    assert field.ctype == q.ElemTypeRealD
    assert isinstance(sf2, q.FieldRealD)
    assert isinstance(mask, q.FieldRealD)
    field._cc_less_than_double(sf2, mask)

def invert_double(field):
    assert field.ctype == q.ElemTypeRealD
    field._cc_invert_double()

def multiply_double(field, factor):
    assert isinstance(field, q.FieldBase)
    assert isinstance(factor, q.FieldBase)
    assert factor.ctype is q.ElemTypeRealD
    field._cc_multiply_double(factor)

def invert_double_field(field):
    """
    cqlat-compatible name for ``invert_double``.
    """
    invert_double(field)

def multiply_double_field(field, factor):
    """
    cqlat-compatible name for ``multiply_double``.
    """
    multiply_double(field, factor)
