import sys
import qlat as q
import numpy as np
from pprint import pformat

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 3],
        [1, 1, 1, 4],
        [1, 1, 1, 6],
        [1, 1, 1, 8],
        [1, 2, 2, 4],
        [2, 2, 2, 4],
        [2, 2, 2, 4],
        ]

q.begin_with_mpi(size_node_list)

if len(sys.argv) == 1:
    q.displayln_info("Usage: topo-measure-wilson-flow [ --source source_config ] [ --output path-for-results ] [ --density-field ] [ --flow-time 6 ] [ --flow-n-step-per-unit-time 80 ]")

q.displayln_info("Topological charge measurement based on Wilson flow with Qlattice")
q.displayln_info("by Luchang Jin")
q.displayln_info("2024/05/20")

p_source = q.get_arg("--source")
p_output = q.get_arg("--output")
p_flow_time = q.get_arg("--flow-time", "6")
p_flow_n_step = q.get_arg("--flow-n-step-per-unit-time", "80")
is_density_field = q.get_option("--density-field")

flow_time = int(p_flow_time)
flow_n_step = int(p_flow_n_step)

info_path = p_output

if is_density_field:
    assert p_output is not None
    density_field_path=info_path
else:
    density_field_path=None

smear_info_list = [
        [ 1.0 / flow_n_step, flow_n_step, 0.0, "runge-kutta", ],
        ] * flow_time

energy_derivative_info = [ 1.0 / flow_n_step, 0.0, "runge-kutta", ]

def load():
    if p_source is None:
        q.displayln_info("Need to provide source file with '--source filename'. Use a sample gauge field for now.")
        total_site = q.Coordinate([ 4, 4, 4, 8, ])
        geo = q.Geometry(total_site)
        gf = q.GaugeField(geo)
        rs = q.RngState("seed")
        gf.set_rand(rs.split("gf-init"), 0.5, 10)
    else:
        gf = q.GaugeField()
        gf.load(p_source)
    q.clear_all_caches()
    q.clear_mem_cache()
    return gf

gf = load()

topo_list, energy_list, = q.smear_measure_topo(
        gf,
        smear_info_list=smear_info_list,
        energy_derivative_info=energy_derivative_info,
        info_path=info_path,
        density_field_path=density_field_path,
        )

if info_path is None:
    q.displayln_info("To save the result, use '--output path'. Print to screen for now.")
    q.displayln_info(pformat(topo_list))
    q.displayln_info(pformat(energy_list))

q.timer_display()

q.end_with_mpi()

exit()
