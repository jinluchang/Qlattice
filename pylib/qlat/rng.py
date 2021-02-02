import cqlat as c

class RngState:

    def __init__(self, v1=None, v2=None):
        if v1 == None:
            self.cdata = c.mk_rng()
        elif type(v1) == RngState and v2 != None:
            rng = v1
            seed = str(v2)
            self.cdata = c.mk_rng(rng.cdata, seed)
        elif type(v1) != RngState and v2 == None:
            seed = str(v1)
            self.cdata = c.mk_rng(rng_state_root.cdata, seed)
        else:
            raise Exception("RngState init")

    def __del__(self):
        c.free_rng(self.cdata)

rng_state_root = RngState()
