#!/usr/bin/env python3

json_results = []
check_eps = 1e-5

import qlat as q
import os
import numpy as np

q.begin_with_mpi()

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)
rs = q.RngState("seed")

n_points = 16
psel = q.PointsSelection()
psel.set_rand(total_site, n_points, rs)

fsel = q.FieldSelection(geo)
fsel.add_psel(psel)

# many checks

q.displayln_info("CHECK: geo.show() =", geo.show())

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))

q.displayln_info(f"CHECK: prop.crc32() = {prop.crc32()} ; prop.qnorm() = {prop.qnorm():.12E}")

prop_arr = np.asarray(prop)
q.displayln_info(f"CHECK: prop_arr.dtype = {prop_arr.dtype}")
q.displayln_info(f"CHECK: prop_arr.shape = {prop_arr.shape}")

prop.save_double("results/prop-double.field")
prop = q.Prop()
prop.load_double("results/prop-double.field")

q.displayln_info(f"CHECK: prop.crc32() = {prop.crc32()} ; prop.qnorm() = {prop.qnorm():.12E}")

prop.save_float_from_double("results/prop-float.field")
prop = q.Prop()
prop.load_double_from_float("results/prop-float.field")

q.displayln_info(f"CHECK: prop.crc32() = {prop.crc32()} ; prop.qnorm() = {prop.qnorm():.12E}")

q.save_pickle_obj(prop, f"results/prop-{q.get_id_node()}.pickle", is_sync_node=False)
prop_load = q.load_pickle_obj(f"results/prop-{q.get_id_node()}.pickle", is_sync_node=False)
assert np.all(prop[:] == prop_load[:])

psel = q.PointsSelection(total_site, [
    [ 0, 0, 0, 0, ],
    [ 0, 1, 2, 0, ],
    ])

q.save_pickle_obj(psel, f"results/psel.pickle")
psel_load = q.load_pickle_obj(f"results/psel.pickle")
assert np.all(psel[:] == psel_load[:])

n_per_tslice = 16
fsel = q.FieldSelection()
fsel.set_rand(geo.total_site, n_per_tslice, rs.split("fsel"))

fselc = fsel.copy()
fselc.add_psel(psel)

q.save_pickle_obj(fselc, f"results/fselc-{q.get_id_node()}.pickle", is_sync_node=False)
fselc_load = q.load_pickle_obj(f"results/fselc-{q.get_id_node()}.pickle", is_sync_node=False)
assert np.all(fselc[:] == fselc_load[:])

s_prop = q.SelProp(fselc)
s_prop @= prop

q.save_pickle_obj(s_prop, f"results/s_prop-{q.get_id_node()}.pickle", is_sync_node=False)
s_prop_load = q.load_pickle_obj(f"results/s_prop-{q.get_id_node()}.pickle", is_sync_node=False)
assert np.all(s_prop[:] == s_prop_load[:])

n1 = s_prop.n_elems
n2 = q.glb_sum(n1)

q.displayln_info(f"CHECK: s_prop n1={n1} ; n2={n2}")

sp_prop = q.PselProp(s_prop, rs.split("shuffle_selected_field-1"))
q.displayln_info(f"CHECK: sp_prop from shuffle_selected_field")

q.save_pickle_obj(sp_prop, f"results/sp_prop-{q.get_id_node()}.pickle", is_sync_node=False)
sp_prop_load = q.load_pickle_obj(f"results/sp_prop-{q.get_id_node()}.pickle", is_sync_node=False)
assert np.all(sp_prop[:] == sp_prop_load[:])

n3 = sp_prop.n_points
n4 = q.glb_sum(n3)
q.displayln_info(f"CHECK: s_prop n3={n3} ; n4={n4}")

assert n4 == n2

q.displayln_info(f"CHECK: len(sp_prop.psel) = {len(sp_prop.psel)}")
q.displayln_info(f"CHECK: sp_prop.psel[:].tolist() = {sp_prop.psel[:].tolist()}")

sig1 = q.get_data_sig(sp_prop[:], q.RngState(f"{q.get_id_node()}"))
sig2 = q.glb_sum(sig1)

json_results.append((f"sp_prop from shuffle_selected_field sig1", sig1, check_eps))
json_results.append((f"sp_prop from shuffle_selected_field sig2", sig2, check_eps))

sp_prop = q.PselProp(psel)
sp_prop @= prop

q.displayln_info("CHECK: sp_prop", sp_prop.qnorm())

sp_prop_arr = np.asarray(sp_prop)
q.displayln_info(f"CHECK: sp_prop_arr.dtype = {sp_prop_arr.dtype}")
q.displayln_info(f"CHECK: sp_prop_arr.shape = {sp_prop_arr.shape}")

sp_prop.save("results/prop.lat")
sp_prop1 = q.PselProp(psel)
sp_prop1.load("results/prop.lat")
sp_prop1 -= sp_prop

q.displayln_info("CHECK: sp_prop", sp_prop.qnorm(), sp_prop1.qnorm(), "save load")

ld = sp_prop.to_lat_data()
sp_prop1 = q.PselProp(psel)
sp_prop1.from_lat_data(ld)
sp_prop1 -= sp_prop

q.displayln_info("CHECK: sp_prop", sp_prop.qnorm(), sp_prop1.qnorm(), "lat_data conversion")

q.save_pickle_obj(sp_prop, f"results/sp_prop.pickle")
sp_prop_load = q.load_pickle_obj(f"results/sp_prop.pickle")
assert np.all(sp_prop[:] == sp_prop_load[:])

s_prop = q.SelProp(fsel)
s_prop @= prop
s_prop.save_double("results/prop.sfield")
s_prop1 = q.SelProp(fsel)
s_prop1.load_double("results/prop.sfield")
s_prop1 -= s_prop

q.displayln_info("CHECK: s_prop", s_prop.qnorm(), s_prop1.qnorm())

s_prop_arr = np.asarray(s_prop)
q.displayln_info(f"CHECK: s_prop_arr.dtype = {s_prop_arr.dtype}")
q.displayln_info(f"CHECK: s_prop_arr.shape = {s_prop_arr.shape}")

prop1 = q.SelProp(fselc)
prop1 @= prop

q.displayln_info(f"CHECK: prop1 {prop1.qnorm():.12E}")

sp_prop1 = q.PselProp(psel)
sp_prop1 @= prop1
sp_prop1 -= sp_prop

q.displayln_info(f"CHECK: sp_prop1 {sp_prop.qnorm():12E} {sp_prop1.qnorm():.12E}")

s_prop1 = q.SelProp(fsel)
s_prop1 @= prop1
s_prop1 -= s_prop

q.displayln_info(f"CHECK: s_prop1 {s_prop.qnorm():12E} {s_prop1.qnorm():.12E}")

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))
q.displayln_info(f"CHECK: prop : {prop.crc32()}")
prop_msc = q.convert_mspincolor_from_wm(prop)
q.displayln_info(f"CHECK: prop_msc : {prop_msc.crc32()}")
prop_wm = q.convert_wm_from_mspincolor(prop_msc)
q.displayln_info(f"CHECK: prop_wm : {prop_wm.crc32()}")
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

q.check_all_files_crc32_info("results")

q.check_log_json(__file__, json_results)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
