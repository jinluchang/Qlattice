import cqlat as c

class GaugeAction:

    def __init__(self, beta, c1 = 0.0):
        self.cdata = c.mk_gauge_action(beta, c1)

    def __del__(self):
        c.free_gauge_action(self)

    def beta(self):
        return c.get_beta_gauge_action(self)

    def c1(self):
        return c.get_c1_gauge_action(self)

