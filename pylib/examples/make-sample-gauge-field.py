#!/usr/bin/env python3

import sys
import qlat as q

q.begin(sys.argv, [
    (1, 1, 1, 1),
    (1, 1, 1, 2),
    (1, 1, 1, 4),
    (1, 1, 1, 8),
    (2, 2, 2, 2),
    (2, 2, 2, 4)])

print("id_node: {:4} / {} ; coor_node: {:9} / {}".format(
    q.get_id_node(),
    q.get_num_node(),
    str(q.get_coor_node()),
    str(q.get_size_node())))

total_site = (4, 4, 4, 8)

geo = q.Geometry(total_site)

if q.get_id_node() == 0:
    print(geo.show_all())

gf = q.GaugeField(geo)

if q.get_id_node() == 0:
    print(gf.geo().show_all())

q.set_unit(gf)

q.gf_show_info(gf)

rng = q.RngState("gf-init")

q.set_g_rand_color_matrix_field(gf, rng, 0.3, 1)

q.gf_show_info(gf)

mview = gf.mview()

if q.get_id_node() == 0:
    print("gf.mview().nbytes =", mview.nbytes)

q.timer_display_stack()

q.timer_display()

q.end()
