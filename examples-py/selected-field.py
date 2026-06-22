#!/usr/bin/env python3

check_eps = 1e-5

import qlat as q
import numpy as np

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.json_results_append(f"selected-field: geo.show()={geo.show()}")
rs = q.RngState("seed")

q.save_pickle_obj(geo, "results/geo.pickle")
geo_load = q.load_pickle_obj("results/geo.pickle")
q.json_results_append(f"selected-field: geo_load.show()={geo_load.show()}")

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))

q.json_results_append(f"prop init crc32 = {prop.crc32()}")
q.json_results_append("prop init qnorm", prop.qnorm(), check_eps)

prop_arr = np.asarray(prop)
q.json_results_append(f"prop_arr.dtype={prop_arr.dtype}")
q.json_results_append(f"prop_arr.shape={prop_arr.shape}")

prop.save_double("results/prop-double.field")
prop = q.Prop()
prop.load_double("results/prop-double.field")

q.json_results_append(f"prop load_double crc32 = {prop.crc32()}")
q.json_results_append("prop load_double qnorm", prop.qnorm(), check_eps)

prop.save_float_from_double("results/prop-float.field")
prop = q.Prop()
prop.load_double_from_float("results/prop-float.field")

q.json_results_append(f"prop load_double_from_float crc32 = {prop.crc32()}")
q.json_results_append("prop load_double_from_float qnorm", prop.qnorm(), check_eps)

q.save_pickle_obj(prop, f"results/prop-{q.get_id_node()}.pickle", is_sync_node=False)
prop_load = q.load_pickle_obj(
    f"results/prop-{q.get_id_node()}.pickle", is_sync_node=False
)
assert np.all(prop[:] == prop_load[:])

psel = q.PointsSelection(
    total_site,
    [
        [
            0,
            0,
            0,
            0,
        ],
        [
            0,
            1,
            2,
            0,
        ],
    ],
)

q.save_pickle_obj(psel, "results/psel.pickle")
psel_load = q.load_pickle_obj("results/psel.pickle")
assert np.all(psel[:] == psel_load[:])

n_per_tslice = 16
fsel = q.FieldSelection()
fsel.set_rand(geo.total_site, n_per_tslice, rs.split("fsel"))

fselc = fsel.copy()
fselc.add_psel(psel)

q.save_pickle_obj(fselc, f"results/fselc-{q.get_id_node()}.pickle", is_sync_node=False)
fselc_load = q.load_pickle_obj(
    f"results/fselc-{q.get_id_node()}.pickle", is_sync_node=False
)
assert np.all(fselc[:] == fselc_load[:])

s_prop = q.SelProp(fselc)
s_prop @= prop

q.save_pickle_obj(
    s_prop, f"results/s_prop-{q.get_id_node()}.pickle", is_sync_node=False
)
s_prop_load = q.load_pickle_obj(
    f"results/s_prop-{q.get_id_node()}.pickle", is_sync_node=False
)
assert np.all(s_prop[:] == s_prop_load[:])

n1 = s_prop.n_elems
n2 = q.glb_sum(n1)

q.json_results_append(f"s_prop n_elems = {n1}")
q.json_results_append(f"s_prop glb_sum n_elems = {n2}")

sp_prop = q.PselProp(
    s_prop,
    q.SelectedShufflePlan(
        "r_from_l",
        q.PointsSelection(s_prop.fsel),
        s_prop.fsel.geo,
        rs.split("shuffle_selected_field-1"),
    ),
)
q.json_results_append("sp_prop from shuffle_selected_field")
q.json_results_append("sp_prop from shuffle_selected_field qnorm", sp_prop.qnorm(), check_eps)

q.save_pickle_obj(
    sp_prop, f"results/sp_prop-{q.get_id_node()}.pickle", is_sync_node=False
)
sp_prop_load = q.load_pickle_obj(
    f"results/sp_prop-{q.get_id_node()}.pickle", is_sync_node=False
)
assert np.all(sp_prop[:] == sp_prop_load[:])

n3 = sp_prop.n_points
n4 = q.glb_sum(n3)
q.json_results_append(f"sp_prop n_points = {n3}")
q.json_results_append(f"sp_prop glb_sum n_points = {n4}")

assert n4 == n2

q.json_results_append("sp_prop psel len", float(len(sp_prop.psel)))
q.json_results_append("sp_prop psel tolist", np.array(sp_prop.psel[:].tolist()))
q.json_results_append("sp_prop psel sig", q.get_data_sig(sp_prop.psel[:], q.RngState()))

sig1 = q.get_data_sig(sp_prop[:], q.RngState(f"{q.get_id_node()}"))
sig2 = q.glb_sum(sig1)

q.json_results_append("sp_prop from shuffle_selected_field sig1", sig1, check_eps)
q.json_results_append("sp_prop from shuffle_selected_field sig2", sig2, check_eps)

sp_prop = q.PselProp(psel)
sp_prop @= prop

q.json_results_append("sp_prop from psel qnorm", sp_prop.qnorm(), check_eps)

sp_prop_arr = np.asarray(sp_prop)
q.json_results_append(f"sp_prop_arr.dtype={sp_prop_arr.dtype}")
q.json_results_append(f"sp_prop_arr.shape={sp_prop_arr.shape}")

sp_prop.save("results/prop.lat")
sp_prop1 = q.PselProp(psel)
sp_prop1.load("results/prop.lat")
sp_prop1 -= sp_prop

q.json_results_append("sp_prop save load qnorm", sp_prop.qnorm(), check_eps)
q.json_results_append("sp_prop save load diff qnorm", sp_prop1.qnorm(), check_eps)

ld = sp_prop.to_lat_data()
sp_prop1 = q.PselProp(psel)
sp_prop1.from_lat_data(ld)
sp_prop1 -= sp_prop

q.json_results_append("sp_prop lat_data qnorm", sp_prop.qnorm(), check_eps)
q.json_results_append("sp_prop lat_data diff qnorm", sp_prop1.qnorm(), check_eps)

q.save_pickle_obj(sp_prop, "results/sp_prop.pickle")
sp_prop_load = q.load_pickle_obj("results/sp_prop.pickle")
assert np.all(sp_prop[:] == sp_prop_load[:])

s_prop = q.SelProp(fsel)
s_prop @= prop
s_prop.save_double("results/prop.sfield")
s_prop1 = q.SelProp(fsel)
s_prop1.load_double("results/prop.sfield")
s_prop1 -= s_prop

q.json_results_append("s_prop save load qnorm", s_prop.qnorm(), check_eps)
q.json_results_append("s_prop save load diff qnorm", s_prop1.qnorm(), check_eps)

s_prop_arr = np.asarray(s_prop)
q.json_results_append(f"s_prop_arr.dtype={s_prop_arr.dtype}")
q.json_results_append(f"s_prop_arr.shape={s_prop_arr.shape}")

prop1 = q.SelProp(fselc)
prop1 @= prop

q.json_results_append("prop1 from s_prop qnorm", prop1.qnorm(), check_eps)

sp_prop1 = q.PselProp(psel)
sp_prop1 @= prop1
sp_prop1 -= sp_prop

q.json_results_append("sp_prop qnorm in sp_prop1 check", sp_prop.qnorm(), check_eps)
q.json_results_append("sp_prop1 diff qnorm", sp_prop1.qnorm(), check_eps)

s_prop1 = q.SelProp(fsel)
s_prop1 @= prop1
s_prop1 -= s_prop

q.json_results_append("s_prop qnorm in s_prop1 check", s_prop.qnorm(), check_eps)
q.json_results_append("s_prop1 diff qnorm", s_prop1.qnorm(), check_eps)

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))
q.json_results_append(f"prop crc32 = {prop.crc32()}")
prop_msc = q.convert_mspincolor_from_wm(prop)
q.json_results_append(f"prop_msc crc32 = {prop_msc.crc32()}")
prop_wm = q.convert_wm_from_mspincolor(prop_msc)
q.json_results_append(f"prop_wm crc32 = {prop_wm.crc32()}")
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

q.check_log_json(__file__)

q.timer_display()

q.end_with_mpi()

q.displayln_info("CHECK: finished successfully.")
