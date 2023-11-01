from qlat_utils import *
from .c import *
from . import c

def set_checkers(field):
    # no longer needed?
    assert field.ctype == ElemTypeRealD
    c.set_checkers_double_field(field)

def set_double_from_complex(field, cf):
    assert isinstance(cf, FieldBase)
    assert field.ctype == ElemTypeRealD
    assert cf.ctype == ElemTypeComplexD
    c.set_double_from_complex_field(field,cf)

def set_complex_from_double(field, sf):
    assert isinstance(sf, FieldBase)
    assert field.ctype == ElemTypeComplexD
    assert sf.ctype == ElemTypeRealD
    c.set_complex_from_double_field(field,sf)

def set_abs_from_complex(field, cf):
    assert isinstance(cf, FieldBase)
    assert field.ctype == ElemTypeRealD
    assert cf.ctype == ElemTypeComplexD
    c.set_abs_from_complex_field(field,cf)

def set_ratio_double(field, sf1, sf2):
    assert isinstance(sf1, FieldBase)
    assert isinstance(sf2, FieldBase)
    assert field.ctype == ElemTypeRealD
    assert sf1.ctype == ElemTypeRealD
    assert sf2.ctype == ElemTypeRealD
    c.set_ratio_double_field(field,sf1,sf2)

def less_than_double(field, sf2, mask):
    assert isinstance(sf2, FieldBase)
    assert isinstance(mask, FieldBase)
    assert field.ctype == ElemTypeRealD
    assert sf2.ctype == ElemTypeRealD
    assert mask.ctype == ElemTypeRealD
    c.less_than_double_field(field,sf2,mask)

def invert_double(field):
    assert field.ctype == ElemTypeRealD
    c.invert_double_field(field)

def multiply_double(field, factor):
    assert isinstance(field, FieldBase)
    assert isinstance(factor, FieldBase)
    assert factor.ctype is ElemTypeRealD
    c.multiply_double_field(field, factor)
