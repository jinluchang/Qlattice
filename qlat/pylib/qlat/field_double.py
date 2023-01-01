from qlat.field import *

def set_checkers(field):
    # no longer needed?
    assert field.ctype == c.ElemTypeDouble
    c.set_checkers_double_field(field)

def set_double_from_complex(field, cf):
    assert isinstance(cf, FieldBase)
    assert field.ctype == c.ElemTypeDouble
    assert cf.ctype == c.ElemTypeComplex
    c.set_double_from_complex_field(field,cf)

def set_complex_from_double(field, sf):
    assert isinstance(sf, FieldBase)
    assert field.ctype == c.ElemTypeComplex
    assert sf.ctype == c.ElemTypeDouble
    c.set_complex_from_double_field(field,sf)

def set_abs_from_complex(field, cf):
    assert isinstance(cf, FieldBase)
    assert field.ctype == c.ElemTypeDouble
    assert cf.ctype == c.ElemTypeComplex
    c.set_abs_from_complex_field(field,cf)

def set_ratio_double(field, sf1, sf2):
    assert isinstance(sf1, FieldBase)
    assert isinstance(sf2, FieldBase)
    assert field.ctype == c.ElemTypeDouble
    assert sf1.ctype == c.ElemTypeDouble
    assert sf2.ctype == c.ElemTypeDouble
    c.set_ratio_double_field(field,sf1,sf2)

def less_than_double(field, sf2, mask):
    assert isinstance(sf2, FieldBase)
    assert isinstance(mask, FieldBase)
    assert field.ctype == c.ElemTypeDouble
    assert sf2.ctype == c.ElemTypeDouble
    assert mask.ctype == c.ElemTypeDouble
    c.less_than_double_field(field,sf2,mask)

def invert_double(field):
    assert field.ctype == c.ElemTypeDouble
    c.invert_double_field(field)
