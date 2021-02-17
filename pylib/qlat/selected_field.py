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
