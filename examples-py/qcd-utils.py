#!/usr/bin/env python3

json_results = []
check_eps = 1e-5

import sys
import qlat as q
import numpy as np

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin_with_mpi(size_node_list)

q.qremove_all_info("results")
q.qmkdir_info("results")

rs = q.RngState("seed")

path = "results/ckpoint.topo1.4nt8.lat"
is_load_config = True
if is_load_config and q.does_file_exist_sync_node(path):
    gf.load(path)
    geo = q.geo_reform(gf.geo)
    total_site = geo.total_site
else:
    total_site = q.Coordinate([ 4, 4, 4, 8, ])
    geo = q.Geometry(total_site)
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split("gf-init"), 0.5, 10)

q.displayln_info(f"CHECK: total_site = {total_site}")
q.displayln_info("CHECK: geo.show() =", geo.show())
gf.show_info()

gf_ape_list = []

gf_ape = gf.copy()
gf_ape_list.append(gf_ape)

gf_ape = q.gf_ape_smear(gf_ape, 0.5, 30)
gf_ape_list.append(gf_ape)

gf_ape = q.gf_spatial_ape_smear(gf_ape, 0.5, 30)
gf_ape_list.append(gf_ape)

for i, gf_ape in enumerate(gf_ape_list):
    plaq = q.gf_avg_plaq(gf_ape)
    spatial_plaq = q.gf_avg_spatial_plaq(gf_ape)
    link_trace = q.gf_avg_link_trace(gf_ape)
    plaq_density = q.gf_plaq_action_density(gf_ape)
    spatial_plaq_density = q.gf_spatial_plaq_action_density(gf_ape)
    topo = q.gf_topology(gf_ape)
    json_results.append((f"gf_ape_list[{i}] plaq", plaq, check_eps,))
    json_results.append((f"gf_ape_list[{i}] plaq2", 1 - plaq_density / 6, check_eps,))
    json_results.append((f"gf_ape_list[{i}] spatial_plaq", spatial_plaq, check_eps,))
    json_results.append((f"gf_ape_list[{i}] spatial_plaq2", 1 - spatial_plaq_density / 3, check_eps,))
    json_results.append((f"gf_ape_list[{i}] link_trace", link_trace, check_eps,))
    json_results.append((f"gf_ape_list[{i}] plaq_density", plaq_density, check_eps,))
    json_results.append((f"gf_ape_list[{i}] spatial_plaq_density", spatial_plaq_density, check_eps,))
    json_results.append((f"gf_ape_list[{i}] topo", topo, check_eps,))

gf_f = gf.copy()

flow_time = 1.0
t = 0

@q.timer_verbose
def measure():
    gf_f.show_info()
    topo_clf = q.gf_topology_clf(gf_f)
    q.displayln_info(f"CHECK: t={t} topo_clf={topo_clf:.10E}")
    topo = q.gf_topology(gf_f)
    topo_terms = q.gf_topology_terms(gf_f)
    topo_field = q.gf_topology_field(gf_f)
    t_sum = topo_field.glb_sum_tslice()
    t_sum = [ str((t, t_sum.get_elem(t).item(),)) for t in range(t_sum.n_points) ]
    q.displayln_info(f"CHECK: t={t} topo_5li={topo:.10E} {sum(topo_terms):.10E}")
    topo_terms_str = ',\n '.join([ str(x) for x in topo_terms ])
    q.displayln_info(f"[ {topo_terms_str},\n]")
    q.displayln_info("\n".join(t_sum))

@q.timer_verbose
def wilson_flow_force(gf, c1 = 0.0):
    ga = q.GaugeAction(3.0, c1)
    gm_force = q.GaugeMomentum()
    q.set_gm_force(gm_force, gf, ga)
    return gm_force

measure()
for i in range(5):
    c1 = -0.331
    q.gf_wilson_flow(gf_f, flow_time, 50, existing_flow_time = t, c1 = c1)
    t += flow_time
    measure()
if False:
    flow_time = 0.2
    for i in range(1000):
        c1 = -1.4008
        force_size = np.sqrt(wilson_flow_force(gf_f, c1).qnorm() / geo.total_volume)
        q.displayln_info(f"CHECK: force_size={force_size:.14E} flow_time={flow_time}")
        q.gf_wilson_flow(gf_f, flow_time, 50, existing_flow_time = t, c1 = c1)
        t += flow_time
        measure()
    gf_f.save("results/ckpoint.topo1.4nt8.lat")

q.smear_measure_topo(gf, is_show_topo_terms=True, density_field_path="results/topo-measure-density")

q.check_log_json(__file__, json_results)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
