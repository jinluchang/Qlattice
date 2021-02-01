import cqlat as c

from qlat.field import *

def set_zero(x):
    if type(x) == Field:
        c.set_zero_field(x.ctype, x.cdata)

def set_unit(x, coef = 1.0):
    if type(x) == Field:
        c.set_unit_field(x.ctype, x.cdata, coef)

