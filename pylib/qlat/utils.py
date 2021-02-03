import cqlat as c

from qlat.field import *

def sqr(x):
    return x * x

def set_zero(x):
    if type(x) == Field:
        c.set_zero_field(x)
    else:
        raise Exception("set_zero")

def set_unit(x, coef = 1.0):
    if type(x) == Field:
        c.set_unit_field(x, coef)
    else:
        raise Exception("set_unit")

def qnorm(x):
    if type(x) == Field:
        c.qnorm_field(x)
    else:
        raise Exception("qnorm")

