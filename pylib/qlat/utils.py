import cqlat as c

from qlat.field import *
from qlat.lat_io import *

from cqlat import index_from_coordinate, coordinate_from_index
from cqlat import random_permute

def sqr(x):
    return x * x

def set_zero(x):
    x.set_zero()

def set_unit(x, coef = 1.0):
    x.set_unit(coef)

def qnorm(x):
    return x.qnorm()

def show(x):
    return x.show()

def unitarize(x):
    x.unitarize()

