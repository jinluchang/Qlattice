import qlat.c as c

from qlat.qcd import *
from qlat.wilson_flow import *

from pprint import pformat

def gf_topology_field_clf(gf : GaugeField):
    geo = gf.geo()
    topf = FieldDouble(geo, 1)
    c.gf_topology_field_clf(topf, gf)
    return topf

@timer
def gf_topology_clf(gf : GaugeField):
    return gf_topology_field_clf(gf).glb_sum()[0]

def gf_topology_field(gf : GaugeField):
    geo = gf.geo()
    topf = FieldDouble(geo, 1)
    c.gf_topology_field(topf, gf)
    return topf

@timer
def gf_topology(gf : GaugeField):
    return gf_topology_field(gf).glb_sum()[0]

def gf_topology_terms_field(gf : GaugeField):
    geo = gf.geo()
    topf = FieldDouble(geo, 5)
    c.gf_topology_terms_field(topf, gf)
    return topf

@timer
def gf_topology_terms(gf : GaugeField):
    return gf_topology_terms_field(gf).glb_sum()

@timer_verbose
def smear_measure_topo(gf, smear_info_list = None, *, is_show_topo_terms = False):
    # smear_info = [ [ step_size, n_step, c1 = 0.0, wilson_flow_integrator_type = "runge-kutta", ], ... ]
    # c1 = 0.0 # Wilson
    # c1 = -0.331 # Iwasaki
    # c1 = -1.4008 # DBW2
    # wilson_flow_integrator_type = "runge-kutta"
    # wilson_flow_integrator_type = "euler"
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
    displayln_info(f"smear_info_list =")
    displayln_info(pformat(smear_info_list))
    flow_time = 0
    topo_list = []
    @timer
    def smear(step_size, n_step, c1 = 0.0, wilson_flow_integrator_type = "runge-kutta"):
        nonlocal flow_time
        flow_time += n_step * step_size
        for i in range(n_step):
            gf_wilson_flow_step(gf, step_size, c1 = c1, wilson_flow_integrator_type = wilson_flow_integrator_type)
    @timer
    def measure():
        gf.show_info()
        plaq = gf.plaq()
        energy_density = gf_energy_density(gf)
        topo_field_clf = gf_topology_field_clf(gf)
        topo_clf = topo_field_clf.glb_sum().item()
        t_sum_clf = topo_field_clf.glb_sum_tslice()
        t_sum_clf = [ t_sum_clf.get_elem(t).item() for t in range(t_sum_clf.n_points()) ]
        topo_field = gf_topology_field(gf)
        topo = topo_field.glb_sum().item()
        t_sum = topo_field.glb_sum_tslice()
        t_sum = [ t_sum.get_elem(t).item() for t in range(t_sum.n_points()) ]
        displayln_info(f"t={flow_time} topo_clf={topo_clf} topo={topo}")
        displayln_info(pformat(list(enumerate(zip(t_sum_clf, t_sum)))))
        topo_list.append({
            "flow_time": flow_time,
            "plaq": plaq,
            "energy_density": energy_density,
            "topo_clf": topo_clf,
            "topo_clf_tslice": t_sum_clf,
            "topo": topo,
            "topo_tslice": t_sum,
            })
        if is_show_topo_terms:
            topo_terms = gf_topology_terms(gf)
            displayln_info(f"t={t} topo={topo} {sum(topo_terms)}")
            topo_terms_str = ',\n '.join([ str(x) for x in topo_terms ])
            displayln_info(f"[ {topo_terms_str},\n]")
    measure()
    for si in smear_info_list:
        smear(*si)
        measure()
    return topo_list
