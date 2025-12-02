# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_types cimport (
    FieldInt,
    FieldRealD,
)
from .geometry cimport Geometry
from .qcd cimport (
    GaugeField,
    GaugeTransform,
)
from .hmc cimport GaugeMomentum
from .gauge_action cimport GaugeAction

import qlat_utils as q
from .hmc import gf_evolve
from .field_utils import (
    shuffle_field,
    shuffle_field_back,
)
from .mpi import (
    glb_sum,
)

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

@q.timer
def gf_block_stout_smear(GaugeField gf, Coordinate block_site, cc.RealD step_size):
    """
    Stout smear in place.
    If block_site = q.Coordinate(), then no blocking.
    Otherwise, perform stout smearing for each block separately.
    """
    cc.gf_block_stout_smear(gf.xxx().val(), gf.xxx().val(), block_site.xx, step_size)

@q.timer
def gf_local_stout_smear(GaugeField gf, Coordinate block_site, cc.RealD step_size):
    """
    Stout smear in place.
    Otherwise, perform stout smearing for local field only.
    No communiation involved.
    If block_site = q.Coordinate(), then no blocking.
    Otherwise, perform stout smearing for each block separately.
    """
    cc.gf_local_stout_smear(gf.xxx().val(), gf.xxx().val(), block_site.xx, step_size)

@q.timer
def gf_local_avg_plaq(GaugeField gf, Coordinate block_site):
    """
    Compuate local (no glb_sum) average of plaq for plaq within with each blocks.
    Used for diagnostic in `gt_local_tree_gauge`.
    If `block_site == q.Coordinate()`, then `block_site = geo.node_site`.
    Otherwise, perform stout smearing for each block separately.
    """
    return cc.gf_local_avg_plaq(gf.xxx().val(), block_site.xx)

@q.timer
def mk_local_tree_gauge_f_dir(
        Geometry geo,
        Coordinate block_site,
        cc.Bool is_uniform,
        RngState rs,
        ):
    """
    return f_dir
    Used for `gt_local_tree_gauge(gf, f_dir)`
    if dir == 5 for a lattice site, then the point is one of the tree root,
    i.e. if f_dir[x] == 5, then gt[x] is identity.
    """
    cdef FieldInt f_dir = FieldInt()
    cc.set_local_tree_gauge_f_dir(f_dir.xx, geo.xx, block_site.xx, is_uniform, rs.xx)
    return f_dir

@q.timer
def gt_local_tree_gauge(GaugeField gf, FieldInt f_dir, cc.Int num_step):
    """
    return gt_inv
    Note that the gauge fixed gauge field should be obtained as
    gf_gt = gt_inv.inv() * gf
    This is opposite to the other gauge tranformation functions.
    is_uniform = True
    block_site = q.Coordinate([ 4, 4, 4, 4, ])
    f_dir = mk_local_tree_gauge_f_dir(gf.geo, block_site, is_uniform, rs)
    """
    cdef Geometry geo = gf.geo
    cdef GaugeTransform gt_inv = GaugeTransform(geo)
    cc.gt_local_tree_gauge(gt_inv.xxx().val(), gf.xxx().val(), f_dir.xx, num_step)
    return gt_inv

@q.timer
def gt_block_tree_gauge(
        GaugeField gf,
        *,
        Coordinate block_site=None,
        Coordinate new_size_node=None,
        cc.Bool is_uniform=True,
        cc.RealD stout_smear_step_size=0.125,
        cc.Int num_smear_step=4,
        list f_dir_list=None,
        RngState rs_f_dir=None,
        ):
    """
    return (gt_inv, f_dir_list,)
    Provide `f_dir_list`, otherwise `rs_f_dir` will be used to generate `f_dir_list`.
    If `f_dir_list` is provided, `rs_f_dir` will be ignored.
    Default `block_site = Coordinate([ 4, 4, 4, 4, ])`
    Default `rs_f_dir = RngState("seed-gt_block_tree_gauge-rs_f_dir")`
    """
    fname = q.get_fname()
    if block_site is None:
        block_site = Coordinate([ 4, 4, 4, 4, ])
    num_step = 0
    for mu in range(4):
        num_step += block_site[mu] // 2
    cdef Geometry geo = gf.geo
    cdef Coordinate total_site = geo.total_site
    if new_size_node is None:
        new_size_node = total_site // block_site
    cdef Coordinate new_node_site = total_site // new_size_node
    assert new_node_site * new_size_node == total_site
    assert (new_node_site // block_site) * block_site == new_node_site
    cdef Coordinate size_node = geo.size_node
    cdef cc.Bool is_shuffle = size_node != new_size_node
    cdef list gf_list
    if is_shuffle:
        gf_list = shuffle_field(gf, new_size_node)
    else:
        gf_list = [ gf.copy(), ]
    cdef list geo_list = [ gf_local.geo for gf_local in gf_list ]
    if f_dir_list is None:
        if rs_f_dir is None:
            rs_f_dir = RngState("seed-gt_block_tree_gauge-rs_f_dir")
        f_dir_list = [
            mk_local_tree_gauge_f_dir(geo_local, block_site, is_uniform, rs_f_dir)
            for geo_local in geo_list
            ]
    cdef list gt_inv_list = []
    cdef GaugeTransform gt_inv
    avg_plaq_sum = 0.0
    avg_plaq_count = 0
    for gf_local, f_dir_local in zip(gf_list, f_dir_list):
        for step in range(num_smear_step):
            gf_local_stout_smear(gf_local, block_site, stout_smear_step_size)
        avg_plaq_count += 1
        avg_plaq_sum = gf_local_avg_plaq(gf_local, block_site)
        gt_inv = gt_local_tree_gauge(gf_local, f_dir_local, num_step)
        gt_inv_list.append(gt_inv)
    gt_inv = GaugeTransform(geo)
    if is_shuffle:
        shuffle_field_back(gt_inv, gt_inv_list, new_size_node)
    else:
        [ gt_inv, ] = gt_inv_list
    avg_plaq = glb_sum(avg_plaq_sum) / glb_sum(avg_plaq_count)
    q.displayln_info(0, f"{fname}: avg_plaq = {avg_plaq}")
    return (gt_inv, f_dir_list,)