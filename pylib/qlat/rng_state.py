import cqlat as c

class RngState:

    def __init__(self, x = None, y = None):
        if x is None:
            self.cdata = c.mk_rng()
        elif isinstance(x, RngState):
            if y is None:
                # make a copy of x
                self.cdata = c.mk_rng(x)
            else:
                seed = str(y)
                self.cdata = c.mk_rng(x, seed)
        else:
            assert y is None
            # seed a new rng
            seed = str(x)
            self.cdata = c.mk_rng(rng_state_root, seed)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_rng(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, RngState)
        c.set_rng(self, v1)
        return self

    def copy(self, is_copying_data = True):
        x = RngState()
        if is_copying_data:
            x @= self
        return x

    def split(self, seed):
        rng = RngState(self, seed)
        return rng

    def rand_gen(self):
        return c.rand_gen(self)

    def u_rand_gen(self, upper = 1.0, lower = 0.0):
        return c.u_rand_gen(self, upper, lower)

    def g_rand_gen(self, center = 0.0, sigma = 1.0):
        return c.g_rand_gen(self, center, sigma)

    def c_rand_gen(self, size):
        # size can be total_site of the lattice
        return c.c_rand_gen(self, size)

    def select(self, l):
        ri = self.rand_gen() % len(l)
        return l[ri]

rng_state_root = RngState()
