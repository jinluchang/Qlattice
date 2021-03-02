#!/usr/bin/env python3

import qlat as q
import os

q.begin()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
geo = q.Geometry((4, 4, 4, 8), 1)
q.displayln_info("geo.show() =", geo.show())

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))

q.displayln_info(f"prop.crc32() = {prop.crc32()} ; prop.qnorm() = {prop.qnorm()}")

prop.save("results/prop.field")
prop = q.Prop()
prop.load("results/prop.field")

q.displayln_info(f"prop.crc32() = {prop.crc32()} ; prop.qnorm() = {prop.qnorm()}")

prop.save_double("results/prop-double.field")
prop = q.Prop()
prop.load_double("results/prop-double.field")

q.displayln_info(f"prop.crc32() = {prop.crc32()} ; prop.qnorm() = {prop.qnorm()}")

prop.save_float_from_double("results/prop-float.field")
prop = q.Prop()
prop.load_double_from_float("results/prop-float.field")

q.displayln_info(f"prop.crc32() = {prop.crc32()} ; prop.qnorm() = {prop.qnorm()}")

psel = q.PointSelection([(0,0,0,0), (0,1,2,0)])
fsel = q.FieldSelection(geo.total_site(), 16, rs.split("fsel"))

fselc = fsel.copy()
fselc.add_psel(psel)

sp_prop = prop.sparse(psel)
sp_prop.save("results/prop.lat")
sp_prop1 = q.PselProp(psel)
sp_prop1.load("results/prop.lat")
sp_prop1 -= sp_prop

q.displayln_info("prop.sparse(psel)", sp_prop.qnorm(), sp_prop1.qnorm())

s_prop = prop.sparse(fsel)
s_prop.save_double("results/prop.sfield")
s_prop1 = q.SelProp(fsel)
s_prop1.load_double("results/prop.sfield")
s_prop1 -= s_prop

q.displayln_info("prop.sparse(fsel)", s_prop.qnorm(), s_prop1.qnorm())

prop1 = prop.sparse(fselc)

q.displayln_info("prop.sparse(fselc)", prop.qnorm())

sp_prop1 = prop1.sparse(psel)
sp_prop1 -= sp_prop

q.displayln_info("prop1.sparse(psel)", sp_prop.qnorm(), sp_prop1.qnorm())

s_prop1 = prop1.sparse(fsel)
s_prop1 -= s_prop

q.displayln_info("prop1.sparse(fsel)", s_prop.qnorm(), s_prop1.qnorm())

if q.get_id_node() == 0:
    q.displayln(os.listdir("results"))

q.timer_display()

q.end()
