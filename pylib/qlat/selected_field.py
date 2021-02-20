import cqlat as c

class SelectedField:

    def __init__(self, ctype, v1 = None, v2 = None, v3 = None):
        assert isinstance(ctype, str)
        self.ctype = ctype
        if v1 == None:
            self.cdata = c.mk_sfield(ctype)
        elif isinstance(v1, Geometry):
            geo = v1
            n_elems = v2
            multiplicity = v3
            assert isinstance(n_elems, int)
            if multiplicity == None:
                self.cdata = c.mk_sfield(ctype, geo, n_elems)
            else:
                assert isinstance(multiplicity, int)
                self.cdata = c.mk_sfield(ctype, geo, n_elems, multiplicity)
        elif isinstance(v1, FieldSelection):
            fsel = v1
            multiplicity = v2
            assert isinstance(multiplicity, int)
            assert v3 == None
            self.cdata = c.mk_sfield_fsel(ctype, fsel, multiplicity)
        else:
            raise Exception("SelectedField.__init__")

    def __del__(self):
        c.free_sfield(self)

    def __imatmul__(self, f1):
        assert isinstance(f1, SelectedField) and f1.ctype == self.ctype
        c.set_sfield(self, f1)
        return self

    def n_elems(self):
        return c.get_n_elems_sfield(self)

    def total_site(self):
        return c.get_total_site_field(self)

    def multiplicity(self):
        return c.get_multiplicity_sfield(self)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_sfield(geo, self)
        return geo

    def __iadd__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_add_sfield(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_sub_sfield(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_sfield(self, factor)
        return self

    def set_zero(self):
        c.set_zero_sfield(self)

    def qnorm(self):
        return c.qnorm_sfield(self)
