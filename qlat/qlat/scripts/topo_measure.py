import sys
import qlat as q
import numpy as np
from pprint import pformat

q.displayln_info("Topological charge measurement with Qlattice")
q.displayln_info("by Luchang Jin")
q.displayln_info("2024/01/25")

if len(sys.argv) == 1:
    q.displayln_info("Usage: topo-measure [ --source source_config ] [ --output output.pickle ] [ --show-topo-terms ] [ --density-field-path path_for_density_field ]")

q.begin_with_mpi()

p_source = q.get_arg("--source")
p_output = q.get_arg("--output")
p_show_topo_terms = q.get_arg("--show-topo-terms")
density_field_path = q.get_arg("--density-field-path")

is_show_topo_terms = p_show_topo_terms is not None

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
    return gf

gf = load()

topo_list = q.smear_measure_topo(gf, is_show_topo_terms=is_show_topo_terms, density_field_path=density_field_path)

if p_output is not None:
    q.save_pickle_obj(topo_list, p_output)
else:
    q.displayln_info("To save the result, use '--output filename.pickle'. Print to screen for now.")
    q.displayln_info(pformat(topo_list))

q.timer_display()

q.end_with_mpi()

exit()
