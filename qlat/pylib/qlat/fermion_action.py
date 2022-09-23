import qlat.cqlat as c

class FermionAction:

    def __init__(self, *, mass, ls, m5, mobius_scale = 1.0, omega = None):
        assert isinstance(mass, float)
        assert isinstance(ls, int)
        assert isinstance(m5, float)
        if omega is None:
            self.cdata = c.mk_fermion_action_mobius(mass, ls, m5, mobius_scale)
        else:
            assert isinstance(omega, list)
            assert ls == len(omega)
            self.cdata = c.mk_fermion_action_zmobius(mass, m5, omega)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_fermion_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, FermionAction)
        c.set_fermion_action(self, v1)
        return self

    def mass(self):
        return c.get_mass_fermion_action(self)

    def ls(self):
        return c.get_ls_fermion_action(self)

    def m5(self):
        return c.get_m5_fermion_action(self)

    def omega(self):
        return c.get_omega_fermion_action(self)

    def mobius_scale(self):
        return c.get_mobius_scale_fermion_action(self)
