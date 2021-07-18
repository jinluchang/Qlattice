import cqlat as c

class RngState:

    def __init__(self, x = None):
        if x == None:
            self.cdata = c.mk_rng()
        elif isinstance(x, RngState):
            # make a copy of x
            self.cdata = c.mk_rng(x)
        elif isinstance(x, str):
            # seed a new rng
            seed = x
            self.cdata = c.mk_rng(rng_state_root, seed)

    def __del__(self):
        c.free_rng(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, RngState)
        c.set_rng(self, v1)
        return self

    def copy(self):
        x = RngState()
        x @= self
        return x

    def split(self, seed):
        rng = RngState()
        rng.cdata = c.mk_rng(self, str(seed))
        return rng

    def rand_gen(self):
        return c.rand_gen(self)

    def u_rand_gen(self, upper = 1.0, lower = 0.0):
        return c.u_rand_gen(self, upper, lower)

    def g_rand_gen(self, center = 0.0, sigma = 1.0):
        return c.g_rand_gen(self, center, sigma)

    def select(self, l):
        ri = self.rand_gen() % len(l)
        return l[ri]

rng_state_root = RngState()
