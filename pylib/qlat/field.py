import cqlat as c

from qlat.geo import *

class Field:

    def __init__(self, ctype, geo = None, multiplicity = None):
        self.ctype = ctype
        if geo == None:
            self.cdata = c.mk_field(ctype)
        elif multiplicity == None:
            self.cdata = c.mk_field(ctype, geo)
        else:
            self.cdata = c.mk_field(ctype, geo, multiplicity)

    def __del__(self):
        c.free_field(self)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_field(geo, self)
        return geo

    def mview(self):
        return c.get_mview_field(self)

    def __imatmul__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_field(self, f1)
        return self

    def __iadd__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_add_field(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_sub_field(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_field(self, factor)
        return self


