#!/usr/bin/env python3

import sys
import qlat as q
import numpy as np

import qlat_gpt as qg

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4]]

# q.begin(sys.argv, size_node_list)

qg.begin_with_gpt()

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = [ 4, 4, 4, 8, ]
geo = q.Geometry(total_site, 1)
q.displayln_info("geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.5, 10)
# gf = qg.load_gauge_field("results/ckpoint.topo1.4nt8.lat")

gf.show_info()

gf_f = gf.copy()

flow_time = 1.0
t = 0

def measure():
    gf_f.show_info()
    topo = q.gf_topology(gf_f)
    topo_terms = q.gf_topology_terms(gf_f)
    topo_field = q.gf_topology_field(gf_f)
    t_sum = topo_field.glb_sum_tslice()
    t_sum = [ str((t, t_sum.get_elem(t).item(),)) for t in range(t_sum.n_points()) ]
    q.displayln_info(f"t={t} topo={topo} {sum(topo_terms)} {topo_terms}")
    q.displayln_info("\n".join(t_sum))

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
        force_size = np.sqrt(wilson_flow_force(gf_f, c1).qnorm() / geo.total_volume())
        q.displayln_info(f"force_size={force_size} flow_time={flow_time}")
        q.gf_wilson_flow(gf_f, flow_time, 50, existing_flow_time = t, c1 = c1)
        t += flow_time
        measure()
    qg.save_gauge_field(gf_f, "results/ckpoint.topo1.4nt8.lat")

q.timer_display()

qg.end_with_gpt()
