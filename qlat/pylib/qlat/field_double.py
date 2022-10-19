from qlat.field import *

def set_checkers(field):
    # no longer needed?
    assert field.ctype == "double"
    c.set_checkers_double_field(field)

def set_double_from_complex(field, cf):
    assert isinstance(cf, Field)
    assert field.ctype == "double"
    assert cf.ctype == "Complex"
    c.set_double_from_complex_field(field,cf)

def set_complex_from_double(field, sf):
    assert isinstance(sf, Field)
    assert field.ctype == "Complex"
    assert sf.ctype == "double"
    c.set_complex_from_double_field(field,sf)

def set_abs_from_complex(field, cf):
    assert isinstance(cf, Field)
    assert field.ctype == "double"
    assert cf.ctype == "Complex"
    c.set_abs_from_complex_field(field,cf)

def set_ratio_double(field, sf1, sf2):
    assert isinstance(sf1, Field)
    assert isinstance(sf2, Field)
    assert field.ctype == "double"
    assert sf1.ctype == "double"
    assert sf2.ctype == "double"
    c.set_ratio_double_field(field,sf1,sf2)

def less_than_double(field, sf2, mask):
    assert isinstance(sf2, Field)
    assert isinstance(mask, Field)
    assert field.ctype == "double"
    assert sf2.ctype == "double"
    assert mask.ctype == "double"
    c.less_than_double_field(field,sf2,mask)

def invert_double(field):
    assert field.ctype == "double"
    c.invert_double_field(field)
