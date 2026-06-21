#!/usr/bin/env python3

import qlat as q
import numpy as np

from qlat_scripts.v1 import (
    set_param,
    get_param,
    get_save_path,
    get_load_path,
    get_job_seed,
)
from qlat_scripts.v1.gen_data import (
    mk_psrc_tag,
    make_fsel_from_weight,
    make_psel_from_weight,
    run_get_inverter,
    run_f_weight_uniform,
    run_f_weight_load,
    run_f_rand_01,
    run_fsel_prob,
    run_psel_prob,
    run_fsel_from_fsel_prob,
    run_psel_from_psel_prob,
    run_fsel_prob_sub_sampling,
    run_psel_prob_sub_sampling,
    run_psel_split,
    run_fsel_split,
    calc_hvp_sum_tslice,
)

q.begin_with_mpi()

q.json_results_append("test gen_data documentation examples")

job_tag = "test-4nt8-docs-gen-data"
traj = 1000

set_param(job_tag, "traj_list")([traj])
set_param(job_tag, "total_site")([4, 4, 4, 8])
set_param(job_tag, "field_selection_fsel_rate")(0.1)
set_param(job_tag, "field_selection_psel_rate")(0.01)

q.json_results_append("mk_psrc_tag")

tag = mk_psrc_tag(q.Coordinate([1, 2, 3, 4]), 0, 1)
assert tag == "xg=(1,2,3,4) ; type=0 ; accuracy=1"
q.json_results_append(f"mk_psrc_tag: {tag}")

tag2 = mk_psrc_tag(q.Coordinate([0, 0, 0, 0]), 1, 2)
assert tag2 == "xg=(0,0,0,0) ; type=1 ; accuracy=2"
q.json_results_append(f"mk_psrc_tag: {tag2}")

q.json_results_append("run_f_weight_uniform")

total_site = q.Coordinate(get_param(job_tag, "total_site"))
geo = q.Geometry(total_site)

get_f_weight = run_f_weight_uniform(job_tag, traj)
assert get_f_weight is not None
f_weight = get_f_weight()
assert f_weight.geo.total_site == total_site
avg = np.average(f_weight[:].ravel())
assert abs(avg - 1.0) < 1e-14
q.json_results_append(f"f_weight uniform: avg={avg}")

q.json_results_append("run_f_weight_load")

get_f_weight_loaded = run_f_weight_load(job_tag, traj)
assert get_f_weight_loaded is not None
f_weight_loaded = get_f_weight_loaded()
assert np.allclose(f_weight[:], f_weight_loaded[:])
q.json_results_append(f"f_weight load match: {np.allclose(f_weight[:], f_weight_loaded[:])}")

q.json_results_append("run_f_rand_01")

get_f_rand_01 = run_f_rand_01(job_tag, traj)
assert get_f_rand_01 is not None
f_rand_01 = get_f_rand_01()
assert f_rand_01.geo.total_site == total_site
r_min = np.min(f_rand_01[:].ravel())
r_max = np.max(f_rand_01[:].ravel())
assert r_min >= 0.0
assert r_max < 1.0
q.json_results_append(f"f_rand_01: min={r_min:.6f} max={r_max:.6f}")

q.json_results_append("run_fsel_prob")

get_fsel_prob = run_fsel_prob(
    job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight
)
assert get_fsel_prob is not None
fsel_prob = get_fsel_prob()
fsel = fsel_prob.fsel
n_elems = q.glb_sum(fsel.n_elems)
q.json_results_append(f"fsel n_elems={n_elems}")

q.json_results_append("run_psel_prob")

get_psel_prob = run_psel_prob(
    job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight
)
assert get_psel_prob is not None
psel_prob = get_psel_prob()
psel = psel_prob.psel
q.json_results_append(f"psel n_points={psel.n_points}")

q.json_results_append("run_fsel_from_fsel_prob / run_psel_from_psel_prob")

get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
get_psel = run_psel_from_psel_prob(get_psel_prob)
assert get_fsel is not None
assert get_psel is not None
fsel_extracted = get_fsel()
psel_extracted = get_psel()
assert q.glb_sum(fsel_extracted.n_elems) == n_elems
assert psel_extracted.n_points == psel.n_points
q.json_results_append(f"extracted fsel n_elems={q.glb_sum(fsel_extracted.n_elems)} psel n_points={psel_extracted.n_points}")

q.json_results_append("run_fsel_prob_sub_sampling")

get_fsel_prob_sub = run_fsel_prob_sub_sampling(
    job_tag,
    traj,
    sub_sampling_rate=0.5,
    get_fsel_prob=get_fsel_prob,
    get_f_rand_01=get_f_rand_01,
    get_f_weight=get_f_weight,
)
assert get_fsel_prob_sub is not None
fsel_prob_sub = get_fsel_prob_sub()
n_sub = q.glb_sum(fsel_prob_sub.fsel.n_elems)
q.json_results_append(f"fsel sub n_elems={n_sub} (original={n_elems})")

q.json_results_append("run_psel_prob_sub_sampling")

get_psel_prob_sub = run_psel_prob_sub_sampling(
    job_tag,
    traj,
    sub_sampling_rate=0.5,
    get_psel_prob=get_psel_prob,
    get_f_rand_01=get_f_rand_01,
    get_f_weight=get_f_weight,
)
assert get_psel_prob_sub is not None
psel_prob_sub = get_psel_prob_sub()
q.json_results_append(f"psel sub n_points={psel_prob_sub.psel.n_points} (original={psel.n_points})")

q.json_results_append("run_psel_split")

get_psel_list = run_psel_split(job_tag, traj, get_psel=get_psel, num_piece=2)
assert get_psel_list is not None
psel_list = get_psel_list()
assert len(psel_list) == 2
total_split = sum(p.n_points for p in psel_list)
q.json_results_append(f"psel_split: pieces={len(psel_list)} total_points={total_split}")

q.json_results_append("run_fsel_split")

get_fsel_list = run_fsel_split(job_tag, traj, get_fsel=get_fsel, num_piece=2)
assert get_fsel_list is not None
fsel_list = get_fsel_list()
assert len(fsel_list) == 2
total_fsel_split = sum(f.n_points for f in fsel_list)
q.json_results_append(f"fsel_split: pieces={len(fsel_list)} total_points={total_fsel_split}")

q.json_results_append("calc_hvp_sum_tslice")

chvp_16 = q.FieldComplexD(geo, 16)
chvp_16.set_zero()
chvp_16[:, 0] = 1.0 + 0.5j
ld_hvp_ts = calc_hvp_sum_tslice(chvp_16)
assert ld_hvp_ts.ndim() == 4
assert ld_hvp_ts.dim_name(0) == "t_dir"
assert ld_hvp_ts.dim_name(1) == "t"
assert ld_hvp_ts.dim_name(2) == "mu"
assert ld_hvp_ts.dim_name(3) == "nu"
assert ld_hvp_ts.dim_size(0) == 4
t_size = max(total_site.to_list())
assert ld_hvp_ts.dim_size(1) == t_size
q.json_results_append(f"hvp_ts shape: {ld_hvp_ts.dim_sizes()}")

q.check_log_json(__file__, check_eps=1e-10)

q.timer_display()

q.end_with_mpi()

q.displayln_info("CHECK: finished successfully.")
