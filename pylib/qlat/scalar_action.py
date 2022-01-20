import cqlat as c

class ScalarAction:

    def __init__(self, m_sq, lmbd):
        self.cdata = c.mk_scalar_action(m_sq, lmbd)

    def __del__(self):
        c.free_scalar_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, ScalarAction)
        c.set_scalar_action(self, v1)
        return self

    def m_sq(self):
        return c.get_m_sq_scalar_action(self)

    def lmbd(self):
        return c.get_lmbd_scalar_action(self)

