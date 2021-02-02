import cqlat as c

class RngState:

    def __init__(self, seed = None):
        if seed == None:
            self.cdata = c.mk_rng()
        else:
            self.cdata = c.mk_rng(rng_state_root.cdata, str(seed))

    def __del__(self):
        c.free_rng(self.cdata)

    def split(self, seed):
        rng = RngState()
        rng.cdata = c.mk_rng(self.cdata, str(seed))
        return rng

    def rand_gen(self):
        return c.rand_gen(self.cdata)

    def u_rand_gen(self, upper = 1.0, lower = 0.0):
        return c.u_rand_gen(self.cdata, upper, lower)

    def g_rand_gen(self, center = 0.0, sigma = 1.0):
        return c.g_rand_gen(self.cdata, center, sigma)

rng_state_root = RngState()
