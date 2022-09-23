import qlat.cqlat as c

class GaugeAction:

    def __init__(self, beta, c1 = 0.0):
        self.cdata = c.mk_gauge_action(beta, c1)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_gauge_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, GaugeAction)
        c.set_gauge_action(self, v1)
        return self

    def beta(self):
        return c.get_beta_gauge_action(self)

    def c1(self):
        return c.get_c1_gauge_action(self)

