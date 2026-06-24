#!/usr/bin/env python3

import qlat as q

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.json_results_append(f"fields-io: geo.show()={geo.show()}")

psel = q.PointsSelection(total_site, [[0, 0, 0, 0], [0, 1, 2, 0]])
n_per_tslice = 16
fsel = q.FieldSelection()
fsel.set_rand(geo.total_site, n_per_tslice, rs.split("fsel"))

fselc = fsel.copy()
fselc.add_psel(psel)

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))

q.json_results_append("fields-io: prop qnorm", prop.qnorm(), 1e-10)
q.json_results_append(f"fields-io: prop crc32 = {prop.crc32()}")

s_prop = q.SelProp(fselc)
s_prop @= prop

sfw = q.open_fields(
    "results/prop1.fields",
    "w",
    q.Coordinate(
        [
            1,
            1,
            1,
            1,
        ]
    ),
)

prop.save_float_from_double(sfw, "prop", skip_if_exist=True)

s_prop.save_float_from_double(sfw, "s_prop", skip_if_exist=True)

sfw.close()

prop_1 = q.Prop()
s_prop_1 = q.SelProp(None)
s_prop_2 = q.SelProp(fselc)
s_prop_3 = q.SelProp(fselc)

sfr = q.open_fields("results/prop1.fields", "r")

prop_1.load_double_from_float(sfr, "prop")
s_prop_1.load_double_from_float(sfr, "s_prop")
s_prop_2.load_double_from_float(sfr, "s_prop")
s_prop_3.load_double_from_float(sfr, "s_prop")

sfr.close()

assert q.is_matching_fsel(s_prop.fsel, fselc)
assert q.is_matching_fsel(s_prop_1.fsel, fselc)
assert q.is_matching_fsel(s_prop_2.fsel, fselc)
assert q.is_matching_fsel(s_prop_3.fsel, fselc)

prop_1 -= prop
s_prop_1 -= s_prop
s_prop_2 -= s_prop
s_prop_3 -= s_prop

assert q.qnorm(prop_1) < 1e-10
assert q.qnorm(s_prop_1) < 1e-10
assert q.qnorm(s_prop_2) < 1e-10
assert q.qnorm(s_prop_3) < 1e-10

sfw = q.open_fields(
    "results/prop.fields",
    "w",
    q.Coordinate(
        [
            1,
            1,
            1,
            8,
        ]
    ),
)

sf_list = sorted(q.show_all_shuffled_fields_writer())
q.sync_node()
q.json_results_append("fields-io: show_all_shuffled_fields_writer n", len(sf_list))
for idx, s in enumerate(sf_list):
    q.json_results_append(f"fields-io: show_all_shuffled_fields_writer {idx} {s}")
q.sync_node()

sq_list = sorted(q.show_all_qfile())
q.sync_node()
q.json_results_append("fields-io: show_all_qfile n", len(sq_list))
for idx, s in enumerate(sq_list):
    q.json_results_append(f"fields-io: show_all_qfile {idx} {s}")
q.sync_node()

q.json_results_append(f"fields-io: sfw.new_size_node()={sfw.new_size_node()}")

prop.save_double(sfw, "prop.d")

prop.save_float_from_double(sfw, "prop")

sfw.flush()

s_prop = q.SelProp(fsel)
s_prop @= prop
q.json_results_append("fields-io: s_prop qnorm", s_prop.qnorm(), 1e-10)
s_prop.save_float_from_double(sfw, "s_prop")

prop1 = q.SelProp(fselc)
prop1 @= prop
q.json_results_append("fields-io: prop1 qnorm", prop1.qnorm(), 1e-10)
prop1.save_float_from_double(sfw, "prop1")

sfw.close()

sfr = q.open_fields("results/prop.fields", "r")

sq_list = sorted(q.show_all_qfile())
q.sync_node()
q.json_results_append("fields-io: show_all_qfile after read n", len(sq_list))
for idx, s in enumerate(sq_list):
    q.json_results_append(f"fields-io: show_all_qfile after read {idx} {s}")
q.sync_node()

fns = sfr.list()

q.json_results_append(f"fields-io: sfr.list()={sfr.list()}")

q.json_results_append(f"fields-io: sfr.new_size_node()={sfr.new_size_node()}")

prop_d = q.Prop()
prop_d.load_double(sfr, "prop.d")
q.json_results_append("fields-io: prop_d qnorm", prop_d.qnorm(), 1e-10)
q.json_results_append(f"fields-io: prop_d crc32 = {prop_d.crc32()}")
prop_d -= prop
q.json_results_append("fields-io: prop_d diff qnorm", prop_d.qnorm(), 1e-10)
q.json_results_append(f"fields-io: prop_d diff crc32 = {prop_d.crc32()}")

prop_f = q.Prop()
prop_f.load_double_from_float(sfr, "prop")
q.json_results_append("fields-io: prop_f qnorm", prop_f.qnorm(), 1e-10)
q.json_results_append(f"fields-io: prop_f crc32 = {prop_f.crc32()}")
prop_f -= prop
q.json_results_append("fields-io: prop_f diff qnorm", prop_f.qnorm(), 1e-10)
q.json_results_append(f"fields-io: prop_f diff crc32 = {prop_f.crc32()}")

s_prop_f = q.SelProp(fsel)
s_prop_f.load_double_from_float(sfr, "s_prop")
q.json_results_append("fields-io: s_prop_f qnorm", s_prop_f.qnorm(), 1e-10)
s_prop_f -= s_prop
q.json_results_append("fields-io: s_prop_f diff qnorm", s_prop_f.qnorm(), 1e-10)

prop1_f = q.SelProp(fselc)
prop1_f.load_double_from_float(sfr, "prop1")
q.json_results_append("fields-io: prop1_f qnorm", prop1_f.qnorm(), 1e-10)
prop1_f -= prop1
q.json_results_append("fields-io: prop1_f diff qnorm", prop1_f.qnorm(), 1e-10)

prop1_ff = q.SelProp(None)
prop1_ff.load_double_from_float(sfr, "prop1")
q.json_results_append("fields-io: prop1_ff qnorm", prop1_ff.qnorm(), 1e-10)
assert q.is_matching_fsel(prop1_ff.fsel, fselc)
prop1_ff -= prop1
q.json_results_append("fields-io: prop1_ff diff qnorm", prop1_ff.qnorm(), 1e-10)

sfr.close()

index_content = q.qcat_bytes_sync_node("results/prop.fields/index.qar")

q.qremove_info("results/prop.fields/index.qar")

q.fields_build_index("results/prop.fields")

index_content2 = q.qcat_bytes_sync_node("results/prop.fields/index.qar")

assert index_content == index_content2

crc = q.compute_crc32("results/prop.fields/index.qar")

q.json_results_append(f"fields-io: index.qar crc={crc:08X}")

q.json_results_append("fields-io: test read_as_char and write")

sfr = q.open_fields("results/prop.fields", "r")
tags = sfr.list()
q.json_results_append(f"fields-io: tags={tags}")
sfw = q.open_fields(
    "results/prop-copy.fields",
    "w",
    q.Coordinate(
        [
            1,
            1,
            1,
            8,
        ]
    ),
)
for tag in tags:
    obj = sfr.read_as_char(tag)
    q.json_results_append(
        f"fields-io: read_as_char tag='{tag}' type='{type(obj).__name__}'"
    )
    obj.save_direct(sfw, tag, skip_if_exist=True)
tags_sfw = sfw.list()
q.json_results_append(f"fields-io: tags_sfw={tags_sfw}")
sfw.close()
sfr.close()

fn1_list = q.qls_all_sync_node("results/prop.fields")
fn2_list = q.qls_all_sync_node("results/prop-copy.fields")
for fn1, fn2 in zip(fn1_list, fn2_list):
    is_reg = q.is_regular_file_sync_node(fn1)
    assert is_reg == q.is_regular_file_sync_node(fn2)
    if is_reg:
        crc1 = q.compute_crc32(fn1)
        crc2 = q.compute_crc32(fn2)
        q.json_results_append(f"fields-io: crc check '{fn1}' crc={crc1:08X}")
        q.json_results_append(f"fields-io: crc check '{fn2}' crc={crc2:08X}")
        assert crc1 == crc2

q.json_results_append(f"fields-io: list_fields={q.list_fields('results/prop.fields')}")

q.properly_truncate_fields("results/prop.fields")

q.json_results_append(f"fields-io: list_fields={q.list_fields('results/prop.fields')}")

q.truncate_fields("results/prop.fields", fns)

q.json_results_append(f"fields-io: list_fields={q.list_fields('results/prop.fields')}")

q.truncate_fields("results/prop.fields", fns[:-1])

q.json_results_append(f"fields-io: list_fields={q.list_fields('results/prop.fields')}")

q.truncate_fields("results/prop.fields", [])

q.json_results_append(f"fields-io: list_fields={q.list_fields('results/prop.fields')}")

q.json_results_append(f"fields-io: qls_all_sync_node={q.qls_all_sync_node('results')}")

q.check_log_json(__file__)
q.timer_display()

q.end_with_mpi()

q.displayln_info("CHECK: finished successfully.")
