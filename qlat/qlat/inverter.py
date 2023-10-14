from qlat_utils import *
import qlat.c as c

from qlat.propagator import *
from qlat.fermion_action import *
from qlat.qcd import *

cache_inv = mk_cache("inv")

class Inverter:

    pass

class InverterDwfFreeField(Inverter):

    # self.mass
    # self.m5
    # self.momtwist
    # self.timer

    def __init__(self, *, mass, m5 = 1.0, momtwist = None, qtimer = TimerNone()):
        if momtwist is None:
            momtwist = [ 0.0, 0.0, 0.0, 0.0, ]
        self.mass = mass
        self.m5 = m5
        self.momtwist = momtwist
        self.timer = qtimer
        assert isinstance(self.mass, float)
        assert isinstance(self.m5, float)
        assert isinstance(self.momtwist, list)
        assert isinstance(self.timer, (Timer, TimerNone,))

    def __mul__(self, prop_src):
        # prop_src: prop or [ prop, ... ]
        if isinstance(prop_src, Prop):
            self.timer.start()
            prop_sol = free_invert(prop_src, self.mass, self.m5, self.momtwist)
            self.timer.stop()
            return prop_sol
        elif isinstance(prop_src, list):
            return [ self * p for p in prop_src ]
        else:
            raise Exception("InverterDwfFreeField")

class InverterDomainWall(Inverter):

	# self.cdata
	# self.timer

    def __init__(self, *, gf, fa, qtimer = TimerNone()):
        self.cdata = c.mk_inverter_domain_wall(gf, fa)
        self.timer = qtimer
        assert isinstance(self.timer, (Timer, TimerNone,))

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_inverter_domain_wall(self)

    def __mul__(self, prop_src):
        # prop_src: prop or [ prop, ... ]
        if isinstance(prop_src, Prop):
            self.timer.start()
            prop_sol = Prop()
            c.invert_inverter_domain_wall(prop_sol, prop_src, self)
            self.timer.stop()
            return prop_sol
        elif isinstance(prop_src, list):
            return [ self * p for p in prop_src ]
        else:
            raise Exception("InverterDomainWall")

    def stop_rsd(self):
        return c.get_stop_rsd_inverter_domain_wall(self)

    def set_stop_rsd(self, stop_rsd):
        return c.set_stop_rsd_inverter_domain_wall(self, stop_rsd)

    def max_num_iter(self):
        return c.get_max_num_iter_inverter_domain_wall(self)

    def set_max_num_iter(self, max_num_iter):
        return c.set_max_num_iter_inverter_domain_wall(self, max_num_iter)

    def max_mixed_precision_cycle(self):
        return c.get_max_mixed_precision_cycle_inverter_domain_wall(self)

    def set_max_mixed_precision_cycle(self, max_mixed_precision_cycle):
        return c.set_max_mixed_precision_cycle_inverter_domain_wall(self, max_mixed_precision_cycle)

class InverterGaugeTransform(Inverter):

    # self.inverter
    # self.gt
    # self.gt_inv
    # self.timer

    def __init__(self, *, inverter, gt,
            qtimer = TimerNone()):
        self.inverter = inverter
        self.gt = gt
        self.timer = qtimer
        assert isinstance(self.inverter, Inverter)
        assert isinstance(self.gt, GaugeTransform)
        assert isinstance(self.timer, (Timer, TimerNone,))
        self.gt_inv = self.gt.inv()

    def __mul__(self, prop_src):
        assert isinstance(prop_src, (Prop, FermionField4d, list,))
        self.timer.start()
        src = self.gt_inv * prop_src
        sol = self.inverter * src
        prop_sol = self.gt * sol
        self.timer.stop()
        return prop_sol
