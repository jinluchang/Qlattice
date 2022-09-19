#!/usr/bin/env python3

import qlat as q
import os

q.begin()

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = [4, 4, 4, 8]
geo = q.Geometry(total_site, 1)
q.displayln_info("geo.show() =", geo.show())
rs = q.RngState("seed")

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

psel = q.PointSelection([
    [ 0, 0, 0, 0, ],
    [ 0, 1, 2, 0, ],
    ],
    geo = geo,
    )
n_per_tslice = 16
fsel = q.FieldSelection(geo.total_site(), n_per_tslice, rs.split("fsel"))

fselc = fsel.copy()
fselc.add_psel(psel)

sp_prop = q.PselProp(psel)
sp_prop @= prop

q.displayln_info("sp_prop", sp_prop.qnorm())

sp_prop.save("results/prop.lat")
sp_prop1 = q.PselProp(psel)
sp_prop1.load("results/prop.lat")
sp_prop1 -= sp_prop

q.displayln_info("sp_prop", sp_prop.qnorm(), sp_prop1.qnorm(), "save load")

ld = sp_prop.to_lat_data()
sp_prop1 = q.PselProp(psel)
sp_prop1.from_lat_data(ld)
sp_prop1 -= sp_prop

q.displayln_info("sp_prop", sp_prop.qnorm(), sp_prop1.qnorm(), "lat_data conversion")

s_prop = q.SelProp(fsel)
s_prop @= prop
s_prop.save_double("results/prop.sfield")
s_prop1 = q.SelProp(fsel)
s_prop1.load_double("results/prop.sfield")
s_prop1 -= s_prop

q.displayln_info("s_prop", s_prop.qnorm(), s_prop1.qnorm())

prop1 = q.SelProp(fselc)
prop1 @= prop

q.displayln_info("prop1", prop1.qnorm())

sp_prop1 = q.PselProp(psel)
sp_prop1 @= prop1
sp_prop1 -= sp_prop

q.displayln_info("sp_prop1", sp_prop.qnorm(), sp_prop1.qnorm())

s_prop1 = q.SelProp(fsel)
s_prop1 @= prop1
s_prop1 -= s_prop

q.displayln_info("s_prop1", s_prop.qnorm(), s_prop1.qnorm())

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))
q.displayln_info(f"prop : {prop.crc32()}")
prop_msc = q.convert_mspincolor_from_wm(prop)
q.displayln_info(f"prop_msc : {prop_msc.crc32()}")
prop_wm = q.convert_wm_from_mspincolor(prop_msc)
q.displayln_info(f"prop_wm : {prop_wm.crc32()}")
prop_wm -= prop
assert prop_wm.qnorm() == 0

s_prop = q.SelProp(fsel)
s_prop_msc = q.SelProp(fsel)
s_prop_wm = q.SelProp(fsel)
s_prop @= prop
s_prop_msc @= prop_msc
s_prop_wm @= prop
s_prop_msc1 = q.convert_mspincolor_from_wm(s_prop)
s_prop_wm1 = q.convert_wm_from_mspincolor(s_prop_msc)
s_prop_msc1 -= s_prop_msc
s_prop_wm1 -= s_prop_wm
assert s_prop_msc1.qnorm() == 0
assert s_prop_wm1.qnorm() == 0

sp_prop = q.PselProp(psel)
sp_prop_msc = q.PselProp(psel)
sp_prop_wm = q.PselProp(psel)
sp_prop @= prop
sp_prop_msc @= prop_msc
sp_prop_wm @= prop
sp_prop_msc1 = q.convert_mspincolor_from_wm(sp_prop)
sp_prop_wm1 = q.convert_wm_from_mspincolor(sp_prop_msc)
sp_prop_msc1 -= sp_prop_msc
sp_prop_wm1 -= sp_prop_wm
assert sp_prop_msc1.qnorm() == 0
assert sp_prop_wm1.qnorm() == 0

if q.get_id_node() == 0:
    q.displayln(os.listdir("results"))

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end()
