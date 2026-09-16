"""
Module ``qlat.inverter_utils``
===============================\n
Pure-Python inverter classes that do not call C++ functions directly.
``Inverter`` (the extension-type base) and ``InverterDomainWall`` live in
``qlat.inverter``; this module holds ``InverterDwfFreeField``,
``InverterGaugeTransform`` and ``EigSystem``.\n
"""

from qlat_utils import (
        CoordinateD,
        Timer,
        TimerNone,
        )

from .qcd import GaugeTransform
from .propagator import (
        Prop,
        FermionField4d,
        free_invert,
        )
from .inverter import Inverter

## -----

class InverterDwfFreeField(Inverter):
    """
    self.mass
    self.m5
    self.momtwist
    self.timer
    """

    def __init__(self, *, mass, m5=1.0, momtwist=None, qtimer=TimerNone()):
        if momtwist is None:
            momtwist = CoordinateD([ 0.0, 0.0, 0.0, 0.0, ])
        self.mass = mass
        self.m5 = m5
        self.momtwist = momtwist
        self.timer = qtimer
        assert isinstance(self.mass, float)
        assert isinstance(self.m5, float)
        assert isinstance(self.momtwist, CoordinateD)
        assert isinstance(self.timer, (Timer, TimerNone,))

    def __mul__(self, prop_src):
        """
        prop_src: prop or [ prop, ... ]
        """
        if isinstance(prop_src, Prop):
            self.timer.start()
            prop_sol = free_invert(prop_src, self.mass, self.m5, self.momtwist)
            self.timer.stop()
            return prop_sol
        elif isinstance(prop_src, list):
            return [self * p for p in prop_src]
        else:
            raise Exception("InverterDwfFreeField")

## -----

class InverterGaugeTransform(Inverter):
    """
    self.inverter
    self.gt
    self.gt_inv
    self.timer
    """

    def __init__(
        self,
        *,
        inverter,
        gt,
        qtimer=TimerNone(),
    ):
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

## -----

class EigSystem:
    pass

## -----
