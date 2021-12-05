import cqlat as c

from qlat.qcd import *

def gf_ape_smear(gf : GaugeField, alpha : float, steps : int = 1):
    return c.gf_ape_smear(gf, gf, alpha, steps)

def gf_spatial_ape_smear(gf : GaugeField, alpha : float, steps : int = 1):
    # used value: alpha = 0.5, steps = 30
    return c.gf_spatial_ape_smear(gf, gf, alpha, steps)

def gf_hyp_smear(gf : GaugeField, alpha1 : float, alpha2 : float, alpha3 : float):
    # values in paper is 0.75 0.6 0.3
    # 10.1103/PhysRevD.64.034504 Eq(4)
    return c.gf_hyp_smear(gf, gf, alpha1, alpha2, alpha3)
