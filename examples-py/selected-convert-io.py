#!/usr/bin/env python3

check_eps = 1e-12

import qlat as q

q.begin_with_mpi()

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)
multiplicity = 3

rs = q.RngState("seed")

rsi = rs.split(f"{geo.id_node}")

n_points = 16
psel = q.PointsSelection()
psel.set_rand(total_site, n_points, rs)

n_per_tslice = 4
fsel = q.FieldSelection(geo)
fsel.set_rand(total_site, n_per_tslice, rs)
fsel.add_psel(psel)

q.displayln_info("CHECK: geo.show()=", geo.show())
q.displayln_info(f"CHECK: psel.points_dist_type={psel.points_dist_type}")
q.displayln_info(f"CHECK: fsel.n_elems={fsel.n_elems}")

f1 = q.FieldComplexD(geo, multiplicity)
f2 = q.FieldComplexF(geo, multiplicity)
f3 = q.FieldRealD(geo, multiplicity)
f4 = q.FieldRealF(geo, multiplicity)

sf1 = q.SelectedFieldComplexD(fsel, multiplicity)
sf2 = q.SelectedFieldComplexF(fsel, multiplicity)
sf3 = q.SelectedFieldRealD(fsel, multiplicity)
sf4 = q.SelectedFieldRealF(fsel, multiplicity)

sp1 = q.SelectedPointsComplexD(psel, multiplicity)
sp2 = q.SelectedPointsComplexF(psel, multiplicity)
sp3 = q.SelectedPointsRealD(psel, multiplicity)
sp4 = q.SelectedPointsRealF(psel, multiplicity)

f1.set_zero()
f2.set_zero()
f3.set_zero()
f4.set_zero()

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

sf1.set_zero()
sf2.set_zero()
sf3.set_zero()
sf4.set_zero()

q.json_results_append(f"q.get_data_sig(sf1, rs)", q.get_data_sig(sf1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf2, rs)", q.get_data_sig(sf2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf3, rs)", q.get_data_sig(sf3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf4, rs)", q.get_data_sig(sf4, rs), check_eps)

sp1.set_zero()
sp2.set_zero()
sp3.set_zero()
sp4.set_zero()

q.json_results_append(f"q.get_data_sig(sp1, rs)", q.get_data_sig(sp1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp2, rs)", q.get_data_sig(sp2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp3, rs)", q.get_data_sig(sp3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp4, rs)", q.get_data_sig(sp4, rs), check_eps)

f1.set_rand(rs)
f2.set_rand(rs)
f3.set_rand(rs)
f4.set_rand(rs)

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

q.json_results_append(f"q.glb_sum(q.get_data_sig(f1[:], rs))", q.glb_sum(q.get_data_sig(f1[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(f2[:], rs))", q.glb_sum(q.get_data_sig(f2[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(f3[:], rs))", q.glb_sum(q.get_data_sig(f3[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(f4[:], rs))", q.glb_sum(q.get_data_sig(f4[:], rsi)), check_eps)

sf1.set_rand(rs)
sf2.set_rand(rs)
sf3.set_rand(rs)
sf4.set_rand(rs)

q.json_results_append(f"q.get_data_sig(sf1, rs)", q.get_data_sig(sf1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf2, rs)", q.get_data_sig(sf2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf3, rs)", q.get_data_sig(sf3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf4, rs)", q.get_data_sig(sf4, rs), check_eps)

q.json_results_append(f"q.glb_sum(q.get_data_sig(sf1[:], rs))", q.glb_sum(q.get_data_sig(sf1[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(sf2[:], rs))", q.glb_sum(q.get_data_sig(sf2[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(sf3[:], rs))", q.glb_sum(q.get_data_sig(sf3[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(sf4[:], rs))", q.glb_sum(q.get_data_sig(sf4[:], rsi)), check_eps)

sp1.set_rand(rs)
sp2.set_rand(rs)
sp3.set_rand(rs)
sp4.set_rand(rs)

q.json_results_append(f"q.get_data_sig(sp1, rs)", q.get_data_sig(sp1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp2, rs)", q.get_data_sig(sp2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp3, rs)", q.get_data_sig(sp3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp4, rs)", q.get_data_sig(sp4, rs), check_eps)

q.json_results_append(f"q.glb_sum(q.get_data_sig(sp1[:], rs))", q.glb_sum(q.get_data_sig(sp1[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(sp2[:], rs))", q.glb_sum(q.get_data_sig(sp2[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(sp3[:], rs))", q.glb_sum(q.get_data_sig(sp3[:], rsi)), check_eps)
q.json_results_append(f"q.glb_sum(q.get_data_sig(sp4[:], rs))", q.glb_sum(q.get_data_sig(sp4[:], rsi)), check_eps)

f1.set_rand_g(rs)
f2.set_rand_g(rs)
f3.set_rand_g(rs)
f4.set_rand_g(rs)

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

sf1.set_rand_g(rs)
sf2.set_rand_g(rs)
sf3.set_rand_g(rs)
sf4.set_rand_g(rs)

q.json_results_append(f"q.get_data_sig(sf1, rs)", q.get_data_sig(sf1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf2, rs)", q.get_data_sig(sf2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf3, rs)", q.get_data_sig(sf3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf4, rs)", q.get_data_sig(sf4, rs), check_eps)

sp1.set_rand_g(rs)
sp2.set_rand_g(rs)
sp3.set_rand_g(rs)
sp4.set_rand_g(rs)

q.json_results_append(f"q.get_data_sig(sp1, rs)", q.get_data_sig(sp1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp2, rs)", q.get_data_sig(sp2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp3, rs)", q.get_data_sig(sp3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp4, rs)", q.get_data_sig(sp4, rs), check_eps)

q.displayln_info(f"CHECK: f1[:].shape={f1[:].shape}")
q.displayln_info(f"CHECK: f1[:].dtype={f1[:].dtype}")

q.displayln_info(f"CHECK: f2[:].shape={f2[:].shape}")
q.displayln_info(f"CHECK: f2[:].dtype={f2[:].dtype}")

q.displayln_info(f"CHECK: f3[:].shape={f3[:].shape}")
q.displayln_info(f"CHECK: f3[:].dtype={f3[:].dtype}")

q.displayln_info(f"CHECK: f4[:].shape={f4[:].shape}")
q.displayln_info(f"CHECK: f4[:].dtype={f4[:].dtype}")

q.displayln_info(f"CHECK: sf1[:].shape={sf1[:].shape}")
q.displayln_info(f"CHECK: sf1[:].dtype={sf1[:].dtype}")

q.displayln_info(f"CHECK: sf2[:].shape={sf2[:].shape}")
q.displayln_info(f"CHECK: sf2[:].dtype={sf2[:].dtype}")

q.displayln_info(f"CHECK: sf3[:].shape={sf3[:].shape}")
q.displayln_info(f"CHECK: sf3[:].dtype={sf3[:].dtype}")

q.displayln_info(f"CHECK: sf4[:].shape={sf4[:].shape}")
q.displayln_info(f"CHECK: sf4[:].dtype={sf4[:].dtype}")

q.displayln_info(f"CHECK: sp1[:].shape={sp1[:].shape}")
q.displayln_info(f"CHECK: sp1[:].dtype={sp1[:].dtype}")

q.displayln_info(f"CHECK: sp2[:].shape={sp2[:].shape}")
q.displayln_info(f"CHECK: sp2[:].dtype={sp2[:].dtype}")

q.displayln_info(f"CHECK: sp3[:].shape={sp3[:].shape}")
q.displayln_info(f"CHECK: sp3[:].dtype={sp3[:].dtype}")

q.displayln_info(f"CHECK: sp4[:].shape={sp4[:].shape}")
q.displayln_info(f"CHECK: sp4[:].dtype={sp4[:].dtype}")

f1.save_double("results/f1.field")
f1r = q.FieldComplexD()
f1r.load_double("results/f1.field")
f1r -= f1
q.json_results_append(f"f1r.qnorm()", f1r.qnorm(), check_eps)

f3.save_double("results/f3.field")
f3r = q.FieldRealD()
f3r.load_double("results/f3.field")
f3r -= f3
q.json_results_append(f"f3r.qnorm()", f3r.qnorm(), check_eps)

sf1.save_double("results/sf1.sfield")
sf1r = q.SelectedFieldComplexD(fsel)
sf1r.load_double("results/sf1.sfield")
sf1r -= sf1
q.json_results_append(f"sf1r.qnorm()", sf1r.qnorm(), check_eps)

sf3.save_double("results/sf3.sfield")
sf3r = q.SelectedFieldRealD(fsel)
sf3r.load_double("results/sf3.sfield")
sf3r -= sf3
q.json_results_append(f"sf3r.qnorm()", sf3r.qnorm(), check_eps)

sp1.save("results/sp1.lat")
sp1r = q.SelectedPointsComplexD(psel)
sp1r.load("results/sp1.lat")
sp1r -= sp1
q.json_results_append(f"sp1r.qnorm()", sp1r.qnorm(), check_eps)

sp2.save("results/sp2.latf")
sp2r = q.SelectedPointsComplexF(psel)
sp2r.load("results/sp2.latf")
sp2r -= sp2
q.json_results_append(f"sp2r.qnorm()", sp2r.qnorm(), check_eps)

sp3.save("results/sp3.lat")
sp3r = q.SelectedPointsRealD(psel)
sp3r.load("results/sp3.lat")
sp3r -= sp3
q.json_results_append(f"sp3r.qnorm()", sp3r.qnorm(), check_eps)

sp4.save("results/sp4.latf")
sp4r = q.SelectedPointsRealF(psel)
sp4r.load("results/sp4.latf")
sp4r -= sp4
q.json_results_append(f"sp4r.qnorm()", sp4r.qnorm(), check_eps)

sf1 @= f1
sf2 @= f2
sf3 @= f3
sf4 @= f4

q.json_results_append(f"q.get_data_sig(sf1, rs)", q.get_data_sig(sf1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf2, rs)", q.get_data_sig(sf2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf3, rs)", q.get_data_sig(sf3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf4, rs)", q.get_data_sig(sf4, rs), check_eps)

sp1 @= f1
sp2 @= f2
sp3 @= f3
sp4 @= f4

q.json_results_append(f"q.get_data_sig(sp1, rs)", q.get_data_sig(sp1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp2, rs)", q.get_data_sig(sp2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp3, rs)", q.get_data_sig(sp3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp4, rs)", q.get_data_sig(sp4, rs), check_eps)

sp1.set_zero()
sp2.set_zero()
sp3.set_zero()
sp4.set_zero()

q.json_results_append(f"q.get_data_sig(sp1, rs)", q.get_data_sig(sp1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp2, rs)", q.get_data_sig(sp2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp3, rs)", q.get_data_sig(sp3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp4, rs)", q.get_data_sig(sp4, rs), check_eps)

sp1 @= sf1
sp2 @= sf2
sp3 @= sf3
sp4 @= sf4

q.json_results_append(f"q.get_data_sig(sp1, rs)", q.get_data_sig(sp1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp2, rs)", q.get_data_sig(sp2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp3, rs)", q.get_data_sig(sp3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sp4, rs)", q.get_data_sig(sp4, rs), check_eps)

# Keep other values in field unchanged

f1 @= sf1
f2 @= sf2
f3 @= sf3
f4 @= sf4

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

# Keep other values in field unchanged

f1 @= sp1
f2 @= sp2
f3 @= sp3
f4 @= sp4

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

# Keep other values in selected field unchanged

sf1 @= sp1
sf2 @= sp2
sf3 @= sp3
sf4 @= sp4

q.json_results_append(f"q.get_data_sig(sf1, rs)", q.get_data_sig(sf1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf2, rs)", q.get_data_sig(sf2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf3, rs)", q.get_data_sig(sf3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf4, rs)", q.get_data_sig(sf4, rs), check_eps)

f1.set_zero()
f2.set_zero()
f3.set_zero()
f4.set_zero()

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

f1 @= sf1
f2 @= sf2
f3 @= sf3
f4 @= sf4

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

f1.set_zero()
f2.set_zero()
f3.set_zero()
f4.set_zero()

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

f1 @= sp1
f2 @= sp2
f3 @= sp3
f4 @= sp4

q.json_results_append(f"q.get_data_sig(f1, rs)", q.get_data_sig(f1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f2, rs)", q.get_data_sig(f2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f3, rs)", q.get_data_sig(f3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(f4, rs)", q.get_data_sig(f4, rs), check_eps)

sf1.set_zero()
sf2.set_zero()
sf3.set_zero()
sf4.set_zero()

q.json_results_append(f"q.get_data_sig(sf1, rs)", q.get_data_sig(sf1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf2, rs)", q.get_data_sig(sf2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf3, rs)", q.get_data_sig(sf3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf4, rs)", q.get_data_sig(sf4, rs), check_eps)

sf1 @= sp1
sf2 @= sp2
sf3 @= sp3
sf4 @= sp4

q.json_results_append(f"q.get_data_sig(sf1, rs)", q.get_data_sig(sf1, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf2, rs)", q.get_data_sig(sf2, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf3, rs)", q.get_data_sig(sf3, rs), check_eps)
q.json_results_append(f"q.get_data_sig(sf4, rs)", q.get_data_sig(sf4, rs), check_eps)

q.check_all_files_crc32_info("results")
q.check_log_json(__file__)
q.timer_display()

q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
