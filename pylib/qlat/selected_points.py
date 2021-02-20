import cqlat as c

class SelectedPoints:

    def __init__(self, ctype, n_points = -1, multiplicity = 1):
        assert isinstance(n_points, int)
        self.ctype = ctype
        if n_points < 0:
            self.cdata = c.mk_spfield(ctype)
        else:
            assert isinstance(multiplicity, int)
            self.cdata = c.mk_spfield(ctype, n_points, multiplicity)

    def __del__(self):
        c.free_spfield(self)

    def __imatmul__(self, f1):
        assert isinstance(f1, SelectedPoints) and f1.ctype == self.ctype
        c.set_spfield(self, f1)
        return self

    def n_points(self):
        return c.get_n_points_spfield(self)

    def multiplicity(self):
        return c.get_multiplicity_spfield(self)

    def __iadd__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_add_spfield(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, Field) and f1.ctype == self.ctype
        c.set_sub_spfield(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_spfield(self, factor)
        return self

    def set_zero(self):
        c.set_zero_spfield(self)

    def qnorm(self):
        return c.qnorm_spfield(self)
