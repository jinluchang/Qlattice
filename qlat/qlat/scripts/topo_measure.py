import sys
import qlat as q
import numpy as np

if len(sys.argv) == 1:
    q.displayln_info("Usage: topo-measure [ --source source_config ] [ --output output.pickle ]")

q.begin_with_mpi()

p_source = q.get_arg("--source")
p_output = q.get_arg("--output")

def load():
    if p_source is None:
        q.displayln_info("Need to provide source file with '--source filename'. Use a sample gauge field for now.")
        total_site = [ 4, 4, 4, 8, ]
        geo = q.Geometry(total_site, 1)
        gf = q.GaugeField(geo)
        rs = q.RngState("seed")
        gf.set_rand(rs.split("gf-init"), 0.5, 10)
    else:
        gf = q.GaugeField()
        gf.load(p_source)
    return gf

gf = load()

topo_list = q.smear_measure_topo(gf)

if p_output is not None:
    q.save_pickle_obj(topo_list, p_output)
else:
    q.displayln_info("To save the result, use '--output filename.pickle'. Print to screen for now.")
    q.displayln_info(q.pformat(topo_list))

q.timer_display()

q.end_with_mpi()

exit()
