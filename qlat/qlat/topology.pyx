# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_types cimport FieldRealD
from .geometry cimport Geometry
from .qcd cimport GaugeField
from .hmc cimport GaugeMomentum
from .gauge_action cimport GaugeAction

from .hmc import set_gm_force, gf_evolve
from .wilson_flow import gf_wilson_flow_step, gf_energy_density_field, gf_energy_derivative_density_field
import qlat_utils as q

from pprint import pformat

@q.timer
def gf_plaq_action_density_field(GaugeField gf):
    r"""
    return paf
    paf.geo.multiplicity == 1
    \sum_P (1 - 1/3 * Re Tr U_P)
    #
    6 plaqs in the above sum (for each site).
    #
    Action = beta * total_volume * 6 * (1 - avg_plaq)
    Action = beta * total_volume * action_density
    Single instanton action = 8 * sqr(PI) / g^2
    beta = 6/g^2
    """
    paf = FieldRealD()
    cc.clf_plaq_action_density_field(paf.xx, gf.xxx().val())
    return paf

@q.timer
def gf_spatial_plaq_action_density_field(GaugeField gf):
    r"""
    return paf
    paf.geo.multiplicity == 1
    \sum_P(spatial only) (1 - 1/3 * Re Tr U_P)
    """
    paf = FieldRealD()
    cc.clf_spatial_plaq_action_density_field(paf.xx, gf.xxx().val())
    return paf

@q.timer
def gf_plaq_action_density(GaugeField gf):
    r"""
    return pa
    ininstance(pa, float)
    pa = gf_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume
    """
    cdef Geometry geo = gf.geo
    cdef cc.Long total_volume = geo.total_volume
    return gf_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume

@q.timer
def gf_spatial_plaq_action_density(GaugeField gf):
    r"""
    return pa
    ininstance(pa, float)
    pa = gf_spatial_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume
    """
    cdef Geometry geo = gf.geo
    cdef cc.Long total_volume = geo.total_volume
    return gf_spatial_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume

@q.timer
def gf_topology_field_clf(GaugeField gf):
    r"""
    return topf
    topf.geo.multiplicity == 1
    Use the basic gf_clover_leaf_field
    NOT using 5 loop improved definition
    """
    topf = FieldRealD()
    cc.clf_topology_field(topf.xx, gf.xxx().val())
    return topf

@q.timer
def gf_topology_clf(GaugeField gf):
    r"""
    return top
    ininstance(top, float)
    Use the basic gf_clover_leaf_field
    NOT using 5 loop improved definition
    """
    return gf_topology_field_clf(gf).glb_sum()[:].item()

@q.timer
def gf_topology_field(GaugeField gf):
    r"""
    return topf
    topf.geo.multiplicity == 1
    Using the 5 loop improved definition
    https://arxiv.org/pdf/hep-lat/9701012v2.pdf
    """
    topf = FieldRealD()
    cc.clf_topology_field_5(topf.xx, gf.xxx().val())
    return topf

@q.timer
def gf_topology(GaugeField gf):
    r"""
    return top
    ininstance(top, float)
    Using the 5 loop improved definition Eq. (2-7)
    https://arxiv.org/pdf/hep-lat/9701012v2.pdf
    """
    return gf_topology_field(gf).glb_sum()[:].item()

@q.timer
def gf_topology_terms_field(GaugeField gf):
    r"""
    return topf;
    topf.geo.multiplicity == 5
    sum of the 5 terms should equal to gf_topology_field
    """
    topf = FieldRealD()
    cc.clf_topology_field_5_terms(topf.xx, gf.xxx().val())
    return topf

@q.timer
def gf_topology_terms(GaugeField gf):
    r"""
    return top_terms;
    top_terms.shape == (5,)
    top_terms.dtype == np.float64
    sum of the 5 terms should equal to gf_topology
    """
    return gf_topology_terms_field(gf).glb_sum()[0, :]

@q.timer_verbose
def smear_measure_topo(
        gf,
        smear_info_list=None,
        *,
        energy_derivative_info=None,
        info_path=None,
        density_field_path=None,
        is_show_topo_terms=False,
        ):
    r"""
    return topo_list, energy_list,
    gf will be modified in-place.
    smear_info_list = [ [ step_size, n_step, c1=0.0, wilson_flow_integrator_type="runge-kutta", ], ... ]
    energy_derivative_info = [ step_size, c1=0.0, wilson_flow_integrator_type="runge-kutta", ]
    c1 = 0.0 # Wilson
    c1 = -0.331 # Iwasaki
    c1 = -1.4008 # DBW2
    wilson_flow_integrator_type = "runge-kutta"
    wilson_flow_integrator_type = "euler"
    `info_path` can be the same as `density_field_path`.
    """
    fname = q.get_fname()
    if smear_info_list is None:
        smear_info_list = [
                [ 0.05, 20, 0.0, "euler", ],
                [ 0.05, 20, 0.0, "euler", ],
                [ 0.05, 20, 0.0, "euler", ],
                [ 0.01, 50, -1.4008, "euler", ],
                [ 0.01, 50, -1.4008, "euler", ],
                [ 0.01, 50, -1.4008, "euler", ],
                [ 0.01, 50, -1.4008, "euler", ],
                ]
    if energy_derivative_info is None:
        energy_derivative_info = [ 0.0125, 0.0, "runge-kutta", ]
    q.displayln_info(0, f"{fname}: smear_info_list =")
    q.displayln_info(0, pformat(smear_info_list))
    geo = gf.geo
    total_volume = geo.total_volume
    total_site = geo.total_site
    spatial_volume = total_volume / total_site[3]
    flow_time = 0
    topo_list = []
    energy_list = []
    @q.timer
    def smear(step_size, n_step, c1=0.0, wilson_flow_integrator_type="runge-kutta"):
        nonlocal flow_time
        for i in range(n_step):
            gf_wilson_flow_step(gf, step_size, c1=c1, wilson_flow_integrator_type=wilson_flow_integrator_type)
            flow_time += step_size
            plaq = gf.plaq()
            energy_density_field = gf_energy_density_field(gf)
            energy_density = energy_density_field.glb_sum()[:].item() / total_volume
            energy_density_tslice = (energy_density_field.glb_sum_tslice()[:].ravel() / spatial_volume).tolist()
            energy_list.append({
                "flow_time": flow_time,
                "plaq": plaq,
                "energy_density": energy_density,
                "energy_density_tslice": energy_density_tslice,
                })
    @q.timer
    def measure():
        gf.show_info()
        plaq = gf.plaq()
        plaq_action_density_field = gf_plaq_action_density_field(gf)
        plaq_action_density = plaq_action_density_field.glb_sum()[:].item() / total_volume
        plaq_action_density_tslice = (plaq_action_density_field.glb_sum_tslice()[:].ravel() / spatial_volume).tolist()
        energy_density_field = gf_energy_density_field(gf)
        energy_density = energy_density_field.glb_sum()[:].item() / total_volume
        energy_density_tslice = (energy_density_field.glb_sum_tslice()[:].ravel() / spatial_volume).tolist()
        topo_field_clf = gf_topology_field_clf(gf)
        topo_clf = topo_field_clf.glb_sum()[:].item()
        topo_clf_tslice = topo_field_clf.glb_sum_tslice()[:].ravel().tolist()
        topo_field = gf_topology_field(gf)
        topo = topo_field.glb_sum()[:].item()
        topo_tslice = topo_field.glb_sum_tslice()[:].ravel().tolist()
        energy_deriv_density_field = gf_energy_derivative_density_field(
                gf,
                epsilon=energy_derivative_info[0],
                c1=energy_derivative_info[1],
                wilson_flow_integrator_type=energy_derivative_info[2],
                )
        energy_deriv = energy_deriv_density_field.glb_sum()[:].item() / total_volume
        energy_deriv_tslice = (energy_deriv_density_field.glb_sum_tslice()[:].ravel() / spatial_volume).tolist()
        q.displayln_info(0, f"{fname}: t={flow_time} ; energy_density={energy_density} ; energy_deriv={energy_deriv} ; topo={topo}")
        q.displayln_info(0, pformat(list(enumerate(zip(energy_density_tslice, energy_deriv_tslice, topo_tslice)))))
        if is_show_topo_terms:
            topo_terms = gf_topology_terms(gf)
            q.displayln_info(0, f"{fname}: t={flow_time} ; topo={topo} ; {sum(topo_terms)}")
            topo_terms_str = ',\n  '.join([ str(x) for x in topo_terms ])
            q.displayln_info(0, f"[ {topo_terms_str},\n]")
        if density_field_path is not None:
            # plaq_action_density_field.save_double(f"{density_field_path}/smear-step-{len(topo_list)}/plaq_action_density.field")
            # topo_field_clf.save_double(f"{density_field_path}/smear-step-{len(topo_list)}/topo_clf.field")
            energy_density_field.save_double(f"{density_field_path}/smear-step-{len(topo_list)}/energy_density.field")
            topo_field.save_double(f"{density_field_path}/smear-step-{len(topo_list)}/topo.field")
            energy_deriv_density_field.save_double(f"{density_field_path}/smear-step-{len(topo_list)}/energy_deriv_density_field.field")
        topo_list.append({
            "flow_time": flow_time,
            "plaq": plaq,
            "energy_density": energy_density,
            "energy_density_tslice": energy_density_tslice,
            "energy_deriv": energy_deriv,
            "energy_deriv_tslice": energy_deriv_tslice,
            "topo": topo,
            "topo_tslice": topo_tslice,
            "topo_clf": topo_clf,
            "topo_clf_tslice": topo_clf_tslice,
            "plaq_action_density": plaq_action_density,
            "plaq_action_density_tslice": plaq_action_density_tslice,
            })
    measure()
    for si in smear_info_list:
        smear(*si)
        measure()
    if info_path is not None:
        q.save_pickle_obj(topo_list, f"{info_path}/info.pickle")
        q.save_pickle_obj(energy_list, f"{info_path}/energy-list.pickle")
    return topo_list, energy_list,
