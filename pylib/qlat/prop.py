import cqlat as c

from qlat.field import *
from qlat.rng import *

class Propagator4d(Field):
    def __init__(self, geo = None):
        Field.__init__(self, "WilsonMatrix", geo, 1)

