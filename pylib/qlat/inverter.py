import cqlat as c

from qlat.propagator import *
from qlat.fermion_action import *
from qlat.qcd import *

class Inverter:

    pass

class InverterDwfFreeField(Inverter):

    def __init__(self, *, mass,
            m5 = 1.0,
            momtwist = [0.0, 0.0, 0.0, 0.0],
            timer = TimerNone()):
        self.mass = mass
        self.m5 = m5
        self.momtwist = momtwist
        self.timer = timer
        assert isinstance(self.mass, float)
        assert isinstance(self.m5, float)
        assert isinstance(self.momtwist, list)
        assert isinstance(self.timer, Timer)

    def __mul__(self, prop_src):
        if isinstance(prop_src, Propagator4d):
            self.timer.start()
            prop_sol = free_invert(prop_src, self.mass, self.m5, self.momtwist)
            self.timer.stop()
            return prop_sol
        elif isinstance(prop_src, list):
            return [ self * p for p in prop_src ]
        else:
            raise Exception("InverterDwfFreeField")

class InverterDomainWall(Inverter):

    def __init__(self, *, gf, fa,
            timer = TimerNone()):
        self.cdata = c.mk_inverter_domain_wall(gf, fa)
        self.timer = timer

    def __del__(self):
        c.free_inverter_domain_wall(self)

    def __mul__(self, prop_src):
        if isinstance(prop_src, Propagator4d):
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

    def __init__(self, *, inverter, gt,
            timer = TimerNone()):
        self.inverter = inverter
        self.gt = gt
        self.timer = timer
        assert isinstance(self.inverter, Inverter)
        assert isinstance(self.gt, GaugeTransform)
        assert isinstance(self.timer, Timer)
        self.gt_inv = self.gt.inv()

    def __mul__(self, prop_src):
        assert isinstance(prop_src, Propagator4d) or isinstance(prop_src, list)
        self.timer.start()
        src = self.gt_inv * prop_src
        sol = self.inverter * src
        prop_sol = self.gt * sol
        self.timer.stop()
        return prop_sol
