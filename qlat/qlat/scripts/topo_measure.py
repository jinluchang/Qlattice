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

q.begin_with_mpi()

if len(sys.argv) == 1:
    q.displayln_info("Usage: topo-measure [ --source source_config ] [ --output path ] [ --density-field ] [ --show-topo-terms ]")

q.displayln_info("Topological charge measurement with Qlattice")
q.displayln_info("by Luchang Jin")
q.displayln_info("2024/01/25")

p_source = q.get_arg("--source")
p_output = q.get_arg("--output")
is_density_field = q.get_option("--density-field")
is_show_topo_terms = q.get_option("--show-topo-terms")

info_path = p_output

if is_density_field:
    assert p_output is not None
    density_field_path=info_path
else:
    density_field_path=None

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
        info_path=info_path,
        density_field_path=density_field_path,
        is_show_topo_terms=is_show_topo_terms,
        )

if info_path is None:
    q.displayln_info("To save the result, use '--output path'. Print to screen for now.")
    q.displayln_info(pformat(topo_list))
    q.displayln_info(pformat(energy_list))

q.timer_display()

q.end_with_mpi()

exit()
