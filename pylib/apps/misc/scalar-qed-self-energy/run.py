#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat_ext as q
import gpt as g
import qlat_gpt as qg
import pprint
import math
import rbc_ukqcd_params as rup
import rbc_ukqcd as ru
import numpy as np

import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [
            "results",
            "/home/frank/application/Public/Muon-GM2-cc/jobs/final-run/all-analysis-data",
            "/home/luchang/application/Public/Muon-GM2-cc/jobs/final-run/all-analysis-data",
            ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer
def mk_pion_prop(total_site, m_pi):
    geo = q.Geometry(total_site, 1)
    grid = qg.mk_grid(geo)
    c = g.complex(grid)
    c[:] = 0
    xg = [ 0, 0, 0, 0, ]
    xl = geo.coordinate_l_from_g(xg)
    c[ [xg,] ] = 1.0
    f = qg.qlat_from_gpt([c])
    f = q.free_scalar_invert_cfield(f, m_pi)
    for t in range(12):
        for i in range(16):
            q.displayln_info("spatial", t, i, f.get_elems([0, 0, i, t,])[0].real)
    return f

r_scaling_factor = 5.0

@q.timer 
def mk_four_point_func_table(total_site, n_dtype):
    info_list = [
            [ "type", n_dtype, ],
            [ "t", 1 + math.ceil(total_site[3] / 2), ],
            [ "r", 1 + math.ceil(1.0 + r_scaling_factor * math.sqrt(3.0) *
                total_site[0] / 2.0), ],
            [ "em", 4, [ "mm", "tt", "ii", "xx", ], ],
            ]
    ld = q.mk_lat_data(info_list)
    return ld

gev_inv_fm = 0.197326979

r_pi_fm_list = [ 0.0, 0.60, 0.65, 0.66, 0.67, 0.70, ]

def interpolate_r_pi_fm(i):
    # best approximate r_pi_fm_list[i]
    i1 = math.floor(i)
    assert i1 >= 0
    i2 = i1 + 1
    if i2 >= len(r_pi_fm_list):
        return r_pi_fm_list[-1]
    elif i1 < 0:
        return r_pi_fm_list[0]
    v1 = r_pi_fm_list[i1]
    v2 = r_pi_fm_list[i2]
    a1 = i2 - i
    a2 = i - i1
    return a1 * v1 + a2 * v2

@q.timer
def mk_four_point_func_table_ff(total_site, m_pi, ainv_gev, ff_tag = ""):
    fn = f"results/table-{total_site[0]}nt{total_site[3]}-{m_pi}-{ainv_gev}-({ff_tag}).lat"
    if q.does_file_exist_sync_node(fn):
        ld = q.load_lat_data(fn)
    else:
        a_fm = gev_inv_fm / ainv_gev
        r_pi_list = [ r_pi_fm / a_fm for r_pi_fm in r_pi_fm_list ]
        n_dtype = len(r_pi_list)
        ld = mk_four_point_func_table(total_site, n_dtype)
        for dtype, r_pi in enumerate(r_pi_list):
            q.displayln_info(f"mk_four_point_func_table_ff: fn='{fn}' dtype={dtype} r_pi={r_pi}")
            f = q.mk_pion_four_point_field(total_site, m_pi, ff_tag, r_pi)
            q.acc_four_point_func_em(ld, f, dtype, r_scaling_factor)
        q.glb_sum(ld)
        q.mk_file_dirs_info(fn)
        ld.save(fn)
    q.partial_sum_r_four_point_func_em(ld)
    return ld

def check_traj(job_tag, traj):
    path = rup.dict_params[job_tag]["data_path_em"]
    fn = f"{path}/results={traj}/four-point-func-em.lat"
    return q.does_file_exist_sync_node(fn)

def find_trajs(job_tag):
    trajs = []
    for traj in range(500, 3000, 10):
        if check_traj(job_tag, traj):
            trajs.append(traj)
    return trajs

rup.dict_params["48I"]["data_path_em"] = get_load_path("lat-four-point-em/48I-0.00078")
rup.dict_params["48I"]["ainv/gev"] = 1.73
rup.dict_params["48I"]["m_pi"] = 0.135 / rup.dict_params["48I"]["ainv/gev"]
rup.dict_params["48I"]["trajs"] = find_trajs("48I")

rup.dict_params["64I"]["data_path_em"] = get_load_path("lat-four-point-em/64I-0.000678")
rup.dict_params["64I"]["ainv/gev"] = 2.359
rup.dict_params["64I"]["m_pi"] = 0.135 / rup.dict_params["64I"]["ainv/gev"]
rup.dict_params["64I"]["trajs"] = find_trajs("64I")

@q.timer
def get_four_point_em(job_tag, traj):
    path = rup.dict_params[job_tag]["data_path_em"]
    fn = f"{path}/results={traj}/four-point-func-em.lat"
    ld = q.LatData()
    ld.load(fn)
    total_site = rup.dict_params[job_tag]["total_site"]
    space_volume = total_site[0] * total_site[1] * total_site[2]
    ld[(0,)] *= 1 / space_volume
    ld[(1,)] *= 1 / space_volume
    return ld

@q.timer
def partial_sum_and_normalize_four_point_em(ld):
    ld = ld.copy()
    q.partial_sum_r_four_point_func_em(ld)
    r_max = ld.dim_size(2) - 1
    for t in range(ld.dim_size(1)):
        fac = -ld[(1, t, r_max, 1,)]
        if not t == 0:
            fac /= 2
        for dtype in [ 0, 1, ]:
            ld[(dtype, t,)] /= fac
    for t in range(ld.dim_size(1)):
        fac = 0
        for dtype in [ 2, 3, 4, 5, ]:
            fac += 0.5 * ld[(dtype, t, r_max, 1,)]
        if not t == 0:
            fac /= 2
        for dtype in [ 2, 3, 4, 5, ]:
            ld[(dtype, t,)] /= fac
    return ld

@q.timer
def get_four_point_em_jk_list(job_tag):
    # with 'partial_sum_and_normalize_four_point_em'
    trajs = rup.dict_params[job_tag]["trajs"]
    eps = 1.0
    jk_list = q.jacknife([ q.Data(get_four_point_em(job_tag, traj)) for traj in trajs ])
    return list(map(partial_sum_and_normalize_four_point_em, map(q.Data.get_val, jk_list)))

@q.timer
def combine_dtypes(ld, pion_type):
    info_list = [ [ "type", 1, ], ] + ld.info()[1:]
    ld_c = q.mk_lat_data(info_list)
    if pion_type == "type1":
        ld_c[(0,)] = ld[(0,)]
    elif pion_type == "type2":
        ld_c[(0,)] = -1 * ld[(1,)]
    elif pion_type == "type3":
        ld_c[(0,)] = 1/2 * (ld[(2,)] + ld[(3,)] + ld[(4,)] + ld[(5,)])
    elif pion_type == "charged-pion":
        ld_c[(0,)] = 4/9 * -1 * ld[(1,)] + 5/9 * 1/2 * (ld[(2,)] + ld[(3,)] + ld[(4,)] + ld[(5,)])
    else:
        raise Exception(f"combine_dtypes pion_type={pion_type}")
    return ld_c

def interpolate_t(ld, dtype, t):
    t1 = math.floor(t)
    t2 = t1 + 1
    t_max = ld.dim_size(1) - 1
    if t2 > t_max:
        return ld[(dtype, t_max,)]
    elif t1 < 0:
        return ld[(dtype, 0,)]
    v1 = ld[(dtype, t1,)]
    v2 = ld[(dtype, t2,)]
    a1 = t2 - t
    a2 = t - t1
    return a1 * v1 + a2 * v2

@q.timer
def get_curve(ld, dtype, t_s, curve_tag):
    # t_s may be L / 2
    v = interpolate_t(ld, dtype, t_s)
    r_len = v.shape[0]
    r_max = r_len - 1
    r_list = list(range(0, r_len, 10)) + [ r_max, ]
    if curve_tag == "tt":
        return np.array([ v[r, 1].real for r in r_list ])
    elif curve_tag == "ii":
        return np.array([ v[r, 2].real for r in r_list ])
    elif curve_tag == "xx":
        return np.array([ v[r, 3].real for r in r_list ])
    else:
        raise Exception("get_curve tag='{tag}'")

@q.timer
def get_curves(ld, t_s, curve_tag):
    dtype_len = ld.dim_size(0)
    return [ get_curve(ld, dtype, t_s, curve_tag) for dtype in range(dtype_len) ]

def interpolate_curve(curve_ff_list, i):
    # best approximate curve_ff_list[i]
    # linear interpolate within the range 0 <= i <= len(curve_ff_list) - 1
    # return boundary value if out of range
    i1 = math.floor(i)
    assert i1 >= 0
    i2 = i1 + 1
    if i2 >= len(curve_ff_list):
        return curve_ff_list[-1]
    elif i1 < 0:
        return curve_ff_list[0]
    v1 = curve_ff_list[i1]
    v2 = curve_ff_list[i2]
    a1 = i2 - i
    a2 = i - i1
    return a1 * v1 + a2 * v2

@q.timer
def match_curve(curve, curve_ff_list, eps = 1e-4, n_divide = 10):
    # return best i so that curve_ff_list[i] approximate curve
    def fcn(i):
        return np.linalg.norm(curve - interpolate_curve(curve_ff_list, i))
    val_min = fcn(0)
    i_min = 0
    def find(i0, i1):
        nonlocal val_min, i_min
        if i1 - i0 < eps:
            return (i0 + i1) / 2
        i_list = [ i0 + n / n_divide * (i1 - i0) for n in range(n_divide + 1) ]
        val_min = fcn(i1)
        n_min = n_divide
        for n in range(len(i_list) - 1):
            val = fcn(i_list[n])
            if val < val_min:
                val_min = val
                n_min = n
        if n_min == 0:
            return find(i_list[0], i_list[1])
        elif n_min == n_divide:
            return find(i_list[-2], i_list[-1])
        else:
            return find(i_list[n_min - 1], i_list[n_min + 1])
    low = 0
    high = len(curve_ff_list) - 1
    return find(low, high)

@q.timer
def diagnose_four_point_em(ld):
    # after 'partial_sum_and_normalize_four_point_em'
    r_max = ld.dim_size(2) - 1
    t_s = math.floor(ld.dim_size(1) / 2)
    if False:
        for dtype in range(ld.dim_size(0)):
            for t in range(ld.dim_size(1)):
                q.displayln_info(f"dtype={dtype} t={t} r={r_max}", "tt", ld[(dtype, t, r_max, 1,)], "ii", ld[(dtype, t, r_max, 2,)])
        for dtype in range(ld.dim_size(0)):
            for r in list(range(0, ld.dim_size(2), 5)) + [ r_max, ]:
                q.displayln_info(f"dtype={dtype} t={t_s} r={r}", "tt", ld[(dtype, t_s, r, 1,)], "ii", ld[(dtype, t_s, r, 2,)], "xx", ld[(dtype, t_s, r, 3,)])
        for dtype in range(ld.dim_size(0)):
            s = 0
            for t in range(ld.dim_size(1)):
                s += q.sqr(t) * ld[(dtype, t, r_max, 2,)]
                tt_val = ld[(dtype, t, r_max, 1,)]
                q.displayln_info(f"dtype={dtype} t={t} \\alpha_\\pi partial sum={s} tt={tt_val}")
    for dtype in range(ld.dim_size(0)):
        q.displayln_info(get_curve(ld, dtype, t_s, "tt"))

def fit_r_pi_fm(job_tag, jk_list_four_point_em, pion_type, ff_tag, curve_tag, t_s_ratio):
    total_site = rup.dict_params[job_tag]["total_site"]
    m_pi = rup.dict_params[job_tag]["m_pi"]
    ainv_gev = rup.dict_params[job_tag]["ainv/gev"]
    jk_list = list(map(lambda ld : combine_dtypes(ld, pion_type), jk_list_four_point_em))
    ld_ff = mk_four_point_func_table_ff(total_site, m_pi, ainv_gev, ff_tag)
    t_s = total_site[0] * t_s_ratio
    curve_ff_list = get_curves(ld_ff, t_s, curve_tag)
    ld = q.jk_avg(jk_list)
    jk_list_r_pi_fm = []
    for ld in jk_list:
        curve = get_curve(ld, 0, t_s, curve_tag)
        i_best = match_curve(curve, curve_ff_list)
        r_pi_fm = interpolate_r_pi_fm(i_best)
        jk_list_r_pi_fm.append(r_pi_fm)
    q.displayln_info(f"r_pi_fm = {q.jk_avg(jk_list_r_pi_fm)} +/- {q.jk_err(jk_list_r_pi_fm)} t_s={t_s} job_tag={job_tag} pion_type={pion_type} ff_tag={ff_tag} curve_tag={curve_tag} t_s_ratio={t_s_ratio}")
    # avg plot
    curve = get_curve(jk_list[0], 0, t_s, curve_tag)
    i_best = match_curve(curve, curve_ff_list)
    r_pi_fm = interpolate_r_pi_fm(i_best)
    curve_ff = interpolate_curve(curve_ff_list, i_best)
    fn = f"results/curve-fits/{job_tag}-{pion_type}-{ff_tag}-{curve_tag}-{t_s_ratio}.txt"
    q.mk_file_dirs_info(fn)
    lines = []
    lines.append(f"# r curve-data curve-form-factor-fit curve-scalar-qed")
    for i, (v1, v2, v3,) in enumerate(zip(curve, curve_ff, curve_ff_list[0])):
        lines.append(f"{i} {v1} {v2} {v3}")
    q.qtouch(fn, "\n".join(lines))

@q.timer
def run_job(job_tag):
    jk_list_four_point_em = get_four_point_em_jk_list(job_tag)
    #
    pion_type_list = [ "type2", "charged-pion", ]
    ff_tag_list = [ "pole", "linear", ]
    curve_tag_list = [ "tt", ]
    t_s_ratio_list = [ 1/2, 1/3, 2/3, ]
    #
    for pion_type in pion_type_list:
        for ff_tag in ff_tag_list:
            for curve_tag in curve_tag_list:
                for t_s_ratio in t_s_ratio_list:
                    fit_r_pi_fm(job_tag, jk_list_four_point_em, pion_type, ff_tag, curve_tag, t_s_ratio)
    #
    q.clean_cache()
    q.timer_display()

qg.begin_with_gpt()

q.check_time_limit()

job_tag_list = [ "48I", "64I", ]

for job_tag in job_tag_list:
    run_job(job_tag)

qg.end_with_gpt()
