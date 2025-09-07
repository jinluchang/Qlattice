# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_types cimport FieldRealD
from .geometry cimport Geometry
from .qcd cimport GaugeField
from .hmc cimport GaugeMomentum
from .gauge_action cimport GaugeAction

from .hmc import gf_evolve
import qlat_utils as q

@q.timer
def gf_energy_density(GaugeField gf):
    return cc.gf_energy_density(gf.xxx().val())

@q.timer
def gf_energy_density_field(GaugeField gf):
    cdef Geometry geo = gf.geo
    cdef FieldRealD fd = FieldRealD(geo, 1)
    cc.gf_energy_density_field(fd.xx, gf.xxx().val())
    return fd

@q.timer
def gf_energy_density_dir_field(GaugeField gf):
    r"""
    Similar to `gf_plaq_field`.
    For smeared gauge field:
    `gf_energy_density_dir_field(gf)[:]` \approx `6 * (1 - gf_plaq_field(gf)[:])`
    `gf_energy_density_dir_field(gf)[:].sum(-1, keepdims=True) == gf_energy_density_field(gf)[:]`
    See
    https://arxiv.org/pdf/1006.4518.pdf Eq. (2.1) (Fig. 1) (approximate Eq. (3.1))
    https://arxiv.org/pdf/1203.4469.pdf
    """
    cdef Geometry geo = gf.geo
    cdef FieldRealD fd = FieldRealD(geo, 6)
    cc.gf_energy_density_dir_field(fd.xx, gf.xxx().val())
    return fd

@q.timer
def gf_plaq_flow_force(GaugeField gf, FieldRealD plaq_factor):
    """
    Compute force with plaq dependent beta factor (relative to standard
    `set_wilson_flow_z`).
    `plaq_factor.multiplicity == 6`.
    Check `gf_plaq_field` for the order of plaq.
    """
    cdef GaugeMomentum gm_force = GaugeMomentum()
    cc.set_plaq_flow_z(gm_force.xxx().val(), gf.xxx().val(), plaq_factor.xx)
    return gm_force

# ----------------------

@q.timer
def gf_wilson_flow_force(GaugeField gf, cc.RealD c1=0.0):
    cdef GaugeMomentum gm_force = GaugeMomentum()
    cc.set_wilson_flow_z(gm_force.xxx().val(), gf.xxx().val(), c1)
    return gm_force

@q.timer
def gf_wilson_flow_step(GaugeField gf, cc.RealD epsilon, *, cc.RealD c1=0.0, str wilson_flow_integrator_type=None):
    """
    Modify `gf` in place.
    default: Runge-Kutta scheme
    http://arxiv.org/abs/1006.4518v3
    """
    if wilson_flow_integrator_type is None:
        wilson_flow_integrator_type = "runge-kutta"
    if wilson_flow_integrator_type == "runge-kutta":
        cc.gf_wilson_flow_step(gf.xxx().val(), epsilon, c1=c1)
    elif wilson_flow_integrator_type == "euler":
        cc.gf_wilson_flow_step_euler(gf.xxx().val(), epsilon, c1=c1)

@q.timer
def gf_energy_derivative_density_field(GaugeField gf, *, cc.RealD epsilon=0.0125, cc.RealD c1=0.0, str wilson_flow_integrator_type=None):
    cdef GaugeField gf1 = gf.copy()
    gf_wilson_flow_step(gf1, epsilon, c1=c1, wilson_flow_integrator_type=wilson_flow_integrator_type)
    cdef FieldRealD fd1 = gf_energy_density_field(gf1)
    gf1 @= gf
    gf_wilson_flow_step(gf1, -epsilon, c1=c1, wilson_flow_integrator_type=wilson_flow_integrator_type)
    cdef FieldRealD fd2 = gf_energy_density_field(gf1)
    fd1 -= fd2
    fd1 *= 1 / (2 * epsilon)
    return fd1

@q.timer
def gf_wilson_flow(GaugeField gf, cc.RealD flow_time, cc.Long steps,
        *, cc.RealD c1=0.0, cc.RealD existing_flow_time=0.0, str wilson_flow_integrator_type=None):
    fname = q.get_fname()
    epsilon = flow_time / steps
    energy_density_list = []
    for i in range(steps):
        gf_wilson_flow_step(gf, epsilon, c1=c1, wilson_flow_integrator_type=wilson_flow_integrator_type)
        t = (i + 1) * epsilon + existing_flow_time
        energy_density = gf_energy_density(gf)
        energy_density_list.append(energy_density)
        q.displayln_info(f"{fname}: t={t} ; E={energy_density} ; t^2 E={t*t*energy_density}")
    return energy_density_list

