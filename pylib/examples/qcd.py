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
    print(q.show_geo_all(geo))

gf = q.Field("ColorMatrix", geo, 4)

if q.get_id_node() == 0:
    print(q.show_geo_all(gf.geo()))

q.set_unit(gf)

q.gf_show_info(gf)

rng = q.RngState("gf-init")

q.set_g_rand_color_matrix_field(gf, rng, 0.3, 1)

q.gf_show_info(gf)

q.timer_display_stack()

q.timer_display()

q.end()
