#!/usr/bin/env python3

import qlat as q
import os

q.begin()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
total_site = [4, 4, 4, 8]
geo = q.Geometry(total_site, 1)
q.displayln_info("geo.show() =", geo.show())

psel = q.PointSelection([[0,0,0,0], [0,1,2,0]])
n_per_tslice = 16
fsel = q.FieldSelection(geo.total_site(), n_per_tslice, rs.split("fsel"))

fselc = fsel.copy()
fselc.add_psel(psel)

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))
q.displayln_info("prop", prop.crc32(), prop.qnorm())

sfw = q.open_fields("results/prop.fields", "w", [1,1,1,8])

q.displayln_info("sfw.new_size_node()", sfw.new_size_node())

prop.save_double(sfw, "prop.d")

prop.save_float_from_double(sfw, "prop")

sfw.flush()

s_prop = prop.sparse(fsel)
q.displayln_info("s_prop = prop.sparse(fsel)", s_prop.qnorm())
s_prop.save_float_from_double(sfw, "s_prop")

prop1 = prop.sparse(fselc)
q.displayln_info("prop1 = prop.sparse(fselc)", prop1.qnorm())
prop1.save_float_from_double(sfw, "prop1")

sfw.close()

sfr = q.open_fields("results/prop.fields", "r")

q.displayln_info("sfr.list()", sfr.list())

q.displayln_info("sfr.new_size_node()", sfr.new_size_node())

prop_d = q.Prop()
prop_d.load_double(sfr, "prop.d")
q.displayln_info("prop_d", prop_d.crc32(), prop_d.qnorm())
prop_d -= prop
q.displayln_info("prop_d -= prop", prop_d.crc32(), prop_d.qnorm())

prop_f = q.Prop()
prop_f.load_double_from_float(sfr, "prop")
q.displayln_info("prop_f", prop_f.crc32(), prop_f.qnorm())
prop_f -= prop
q.displayln_info("prop_f -= prop", prop_f.crc32(), prop_f.qnorm())

s_prop_f = q.SelProp(fsel)
s_prop_f.load_double_from_float(sfr, "s_prop")
q.displayln_info("s_prop_f", s_prop_f.qnorm())
s_prop_f -= s_prop
q.displayln_info("s_prop_f -= s_prop", s_prop_f.qnorm())

prop1_f = q.SelProp(fselc)
prop1_f.load_double_from_float(sfr, "prop1")
q.displayln_info("prop1_f", prop1_f.qnorm())
prop1_f -= prop1
q.displayln_info("prop1_f -= prop1", prop1_f.qnorm())

sfr.close()

q.displayln_info(q.list_fields("results/prop.fields"))

q.properly_truncate_fields_sync_node("results/prop.fields")

if q.get_id_node() == 0:
    q.displayln_info(os.listdir("results"))

q.timer_display()

q.end()
