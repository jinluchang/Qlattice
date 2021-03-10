from qlat.propagator import *
from qlat.qcd import *

class Inverter:

    pass

class InverterDwfFreeField(Inverter):

    def __init__(self, **kwargs):
        self.mass = kwargs["mass"]
        self.m5 = kwargs.get("m5", 1.0)
        self.momtwist = kwargs.get("momtwist", [0.0, 0.0, 0.0, 0.0])
        self.timer = kwargs.get("timer", TimerNone())
        assert isinstance(self.mass, float)
        assert isinstance(self.m5, float)
        assert isinstance(self.momtwist, list)
        assert isinstance(self.timer, Timer)

    def __mul__(self, prop_src):
        assert isinstance(prop_src, Propagator4d)
        self.timer.start()
        prop_sol = free_invert(prop_src, self.mass, self.m5, self.momtwist)
        self.timer.stop()
        return prop_sol

class InverterGaugeTransform(Inverter):

    def __init__(self, **kwargs):
        self.inverter = kwargs["inverter"]
        self.gt = kwargs["gt"]
        self.timer = kwargs.get("timer")
        assert isinstance(self.inverter, Inverter)
        assert isinstance(self.gt, GaugeTransform)
        assert isinstance(self.timer, Timer)
        self.gt_inv = self.gt.inv()

    def __mul__(self, prop_src):
        assert isinstance(prop_src, Propagator4d)
        self.timer.start()
        src = self.gt_inv * prop_src
        sol = self.inverter * src
        prop_sol = self.gt * sol
        self.timer.stop()
        return prop_sol
