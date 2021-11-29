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
import scipy as sp

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

def partial_sum_arr(arr):
    size = arr.shape[0]
    arr = np.copy(arr)
    for i in range(1, size):
        arr[i] += arr[i - 1]
    return arr

@q.timer
def partial_sum_r_four_point_func_em(ld):
    sizes = ld.dim_sizes()
    ld = ld.copy()
    arr = ld[()]
    (dtype_size, t_size, r_size, em_size,) = arr.shape
    for dtype in range(dtype_size):
        for t in range(t_size):
            arr[dtype, t] = partial_sum_arr(arr[dtype, t])
    ld[()] = arr
    return ld

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

r_pi_fm_list = [ 0.0, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, ]

def interpolate(v, i):
    size = len(v)
    i1 = math.floor(i)
    assert i1 >= 0
    i2 = i1 + 1
    if i2 >= size:
        return v[size - 1]
    elif i1 < 0:
        return v[0]
    v1 = v[i1]
    v2 = v[i2]
    a1 = i2 - i
    a2 = i - i1
    return a1 * v1 + a2 * v2

def interpolate_r_pi_fm(i):
    # best approximate r_pi_fm_list[i]
    return interpolate(r_pi_fm_list, i)

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
    return ld

def check_traj(job_tag, traj):
    path = rup.dict_params[job_tag]["data_path_em"]
    fn = f"{path}/results={traj}/four-point-func-em.lat"
    fnw = f"{path}/results={traj}/four-point-func-emw.lat"
    return q.does_file_exist_sync_node(fn) or q.does_file_exist_sync_node(fnw)

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

rup.dict_params["24D"]["data_path_em"] = get_load_path("lat-four-point-emw/24D-0.00107")
rup.dict_params["24D"]["ainv/gev"] = 1.0158
rup.dict_params["24D"]["m_pi"] = 0.13975
rup.dict_params["24D"]["trajs"] = find_trajs("24D")

rup.dict_params["24DH"]["data_path_em"] = get_load_path("lat-four-point-emw/24D-0.0174")
rup.dict_params["24DH"]["ainv/gev"] = 1.0158
rup.dict_params["24DH"]["m_pi"] = 0.3357
rup.dict_params["24DH"]["trajs"] = find_trajs("24DH")

rup.dict_params["32D"]["data_path_em"] = get_load_path("lat-four-point-emw/32D-0.00107")
rup.dict_params["32D"]["ainv/gev"] = 1.0158
rup.dict_params["32D"]["m_pi"] = 0.139474
rup.dict_params["32D"]["trajs"] = find_trajs("32D")

rup.dict_params["32Dfine"]["data_path_em"] = get_load_path("lat-four-point-emw/32Dfine-0.0001")
rup.dict_params["32Dfine"]["ainv/gev"] = 1.378
rup.dict_params["32Dfine"]["m_pi"] = 0.10468
rup.dict_params["32Dfine"]["trajs"] = find_trajs("32Dfine")

@q.timer
def get_four_point_em(job_tag, traj):
    path = rup.dict_params[job_tag]["data_path_em"]
    fn = f"{path}/results={traj}/four-point-func-em.lat"
    if not q.does_file_exist_sync_node(fn):
        fn = f"{path}/results={traj}/four-point-func-emw.lat"
    ld = q.LatData()
    ld.load(fn)
    total_site = rup.dict_params[job_tag]["total_site"]
    space_volume = total_site[0] * total_site[1] * total_site[2]
    ld[(0,)] *= 1 / space_volume
    ld[(1,)] *= 1 / space_volume
    return ld

@q.timer
def normalize_four_point_em(ld):
    ld = ld.copy()
    ld_sum = partial_sum_r_four_point_func_em(ld)
    r_max = ld.dim_size(2) - 1
    n_dtype = ld.dim_size(0)
    assert n_dtype >= 2
    for t in range(ld.dim_size(1)):
        fac = -ld_sum[(1, t, r_max, 1,)]
        if not t == 0:
            fac /= 2
        for dtype in [ 0, 1, ]:
            ld[(dtype, t,)] /= fac
    if n_dtype >= 6:
        for t in range(ld.dim_size(1)):
            fac = 0
            for dtype in [ 2, 3, 4, 5, ]:
                fac += 0.5 * ld_sum[(dtype, t, r_max, 1,)]
            if not t == 0:
                fac /= 2
            for dtype in [ 2, 3, 4, 5, ]:
                ld[(dtype, t,)] /= fac
    return ld

@q.timer
def get_four_point_em_jk_list(job_tag):
    # with 'normalize_four_point_em'
    trajs = rup.dict_params[job_tag]["trajs"]
    eps = 1.0
    jk_list = q.jackknife([ q.Data(get_four_point_em(job_tag, traj)) for traj in trajs ])
    return list(map(normalize_four_point_em, map(q.Data.get_val, jk_list)))

@q.timer
def combine_dtypes(ld, pion_type):
    info_list = [ [ "type", 1, ], ] + ld.info()[1:]
    ld_c = q.mk_lat_data(info_list)
    n_dtype = ld.dim_size(0)
    if pion_type == "type1":
        ld_c[(0,)] = ld[(0,)]
    elif pion_type == "type2":
        ld_c[(0,)] = -1 * ld[(1,)]
    elif pion_type == "type3":
        assert n_dtype >= 6
        ld_c[(0,)] = 1/2 * (ld[(2,)] + ld[(3,)] + ld[(4,)] + ld[(5,)])
    elif pion_type == "charged-pion":
        assert n_dtype >= 6
        ld_c[(0,)] = 4/9 * -1 * ld[(1,)] + 5/9 * 1/2 * (ld[(2,)] + ld[(3,)] + ld[(4,)] + ld[(5,)])
    elif pion_type == "pion-diff":
        ld_c[(0,)] = ld[(0,)] - ld[(1,)]
    else:
        raise Exception(f"combine_dtypes pion_type={pion_type}")
    return ld_c

def interpolate_t(ld, dtype, t):
    arr = ld[(dtype,)]
    return interpolate(arr, t)

def partial_sum_with_r_sqr_flat(arr):
    size = arr.shape[0]
    arr_original = arr
    arr_rsum = np.copy(arr)
    for ri in reversed(range(size - 1)):
        arr_rsum[ri] += arr_rsum[ri + 1]
    arr = np.copy(arr)
    s = arr[0] * 0
    for ri in range(1, size):
        r = ri / r_scaling_factor
        s += arr[ri] * r**2
        if ri + 1 < size:
            arr[ri] = s + arr_rsum[ri + 1] * r**2
        else:
            arr[ri] = s
    return arr

def partial_sum_with_r_sqr(arr):
    size = arr.shape[0]
    arr_original = arr
    arr = np.copy(arr)
    s = arr[0] * 0
    for ri in range(1, size):
        r = ri / r_scaling_factor
        s += arr[ri] * r**2
        arr[ri] = s
    return arr

def partial_sum_with_r4(arr):
    size = arr.shape[0]
    arr_original = arr
    arr = np.copy(arr)
    s = arr[0] * 0
    for ri in range(1, size):
        r = ri / r_scaling_factor
        s += arr[ri] * r**4
        arr[ri] = s
    return arr

@q.timer
def get_curve(ld, dtype, t_s, a_fm, curve_tag):
    r_range_tag, em_tag = curve_tag
    # t_s may be L / 2
    v = interpolate_t(ld, dtype, t_s)
    v_list = []
    # ADJUST ME
    # v_list.append(partial_sum_arr(v))
    v_list.append(partial_sum_with_r_sqr(v)) # default
    v_list.append(partial_sum_with_r4(v)) # default
    # v_list.append(partial_sum_with_r_sqr_flat(v))
    #
    if r_range_tag == "all":
        r_len = v.shape[0]
        r_max = r_len - 1
        r_list = list(range(0, r_len, 5)) + [ r_max, ]
    else:
        r_cut_i = float(r_range_tag) / a_fm * r_scaling_factor
        r_list = [ r_cut_i, ]
    if em_tag == "tt":
        em = 1
    elif em_tag == "ii":
        em = 2
    elif em_tag == "xx":
        em = 3
    else:
        raise Exception("get_curve tag='{tag}'")
    return np.array([ interpolate(v, r)[em].real for v in v_list for r in r_list ])

@q.timer
def get_curves(ld, t_s, a_fm, curve_tag):
    dtype_len = ld.dim_size(0)
    return [ get_curve(ld, dtype, t_s, a_fm, curve_tag) for dtype in range(dtype_len) ]

def interpolate_curves(curve_ff_list, i):
    # best approximate curve_ff_list[i]
    # linear interpolate within the range 0 <= i <= len(curve_ff_list) - 1
    # return boundary value if out of range
    return interpolate(curve_ff_list, i)

@q.timer
def match_curve(curve, curve_ff_list_list, *, eps = 1e-4, n_divide = 10):
    # return best i so that \sum_k a_k curve_ff_list_list[k][i] (where \sum_k a_k = 1) approximate curve
    # curve_ff_list_list = [ curve_ff_list_for_one_ff_tag, ... ]
    # curve_ff_list_for_one_ff_tag = [ curve_for_one_r_pi, ... ]
    # curve_for_one_r_pi = [ curve_value_at_one_r, ... ]
    def fcn(i):
        curve_i_list = [ interpolate_curves(curve_ff_list, i) for curve_ff_list in curve_ff_list_list ]
        curve_i = curve_i_list[0]
        return np.linalg.norm(curve - curve_i)
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
    high = len(curve_ff_list_list[0]) - 1
    return find(low, high)

@q.timer
def diagnose_four_point_em(ld):
    # after 'normalize_four_point_em'
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
    a_fm = 0.1
    for dtype in range(ld.dim_size(0)):
        q.displayln_info(get_curve(ld, dtype, t_s, a_fm, "tt"))

def fit_r_pi_fm(job_tag, jk_list_four_point_em, pion_type, ff_tag_list, curve_tag, t_s_fm):
    total_site = rup.dict_params[job_tag]["total_site"]
    m_pi = rup.dict_params[job_tag]["m_pi"]
    ainv_gev = rup.dict_params[job_tag]["ainv/gev"]
    a_fm = gev_inv_fm / ainv_gev
    jk_list = list(map(lambda ld : combine_dtypes(ld, pion_type), jk_list_four_point_em))
    t_s = t_s_fm / a_fm
    curve_ff_list_list = []
    for ff_tag in ff_tag_list:
        ld_ff = mk_four_point_func_table_ff(total_site, m_pi, ainv_gev, ff_tag)
        curve_ff_list_list.append(get_curves(ld_ff, t_s, a_fm, curve_tag))
    ld = q.jk_avg(jk_list)
    jk_list_r_pi_fm = []
    for ld in jk_list:
        curve = get_curve(ld, 0, t_s, a_fm, curve_tag)
        i_best = match_curve(curve, curve_ff_list_list)
        r_pi_fm = interpolate_r_pi_fm(i_best)
        jk_list_r_pi_fm.append(r_pi_fm)
    result_str = f"r_pi_fm = {q.jk_avg(jk_list_r_pi_fm)} +/- {q.jk_err(jk_list_r_pi_fm)} t_s={t_s} job_tag={job_tag} pion_type={pion_type} ff_tag_list={ff_tag_list} curve_tag={curve_tag} t_s_fm={t_s_fm}"
    q.displayln_info(result_str)
    # avg plot
    curve = get_curve(jk_list[0], 0, t_s, a_fm, curve_tag)
    i_best = match_curve(curve, curve_ff_list_list)
    r_pi_fm = interpolate_r_pi_fm(i_best)
    curve_ff = interpolate_curves(curve_ff_list_list[0], i_best)
    fn = f"results/curve-fits/{job_tag}-{pion_type}-{ff_tag_list}-{curve_tag}-{t_s_fm}.txt"
    q.mk_file_dirs_info(fn)
    lines = []
    lines.append(f"# r curve-data curve-form-factor-fit({ff_tag_list[0]}) curve-scalar-qed")
    lines.append(f"# {result_str}")
    for i, (v1, v2, v3,) in enumerate(zip(curve, curve_ff, curve_ff_list_list[0][0])):
        lines.append(f"{i} {v1} {v2} {v3}")
    lines.append("")
    q.qtouch(fn, "\n".join(lines))
    return jk_list_r_pi_fm

@q.timer
def run_job(job_tag):
    jk_list_four_point_em = get_four_point_em_jk_list(job_tag)
    # ADJUST ME
    pion_type_list = [
            "type2",
            "pion-diff",
            # "charged-pion",
            ]
    ff_tag_list_list = [
            [ "pole", ],
            [ "pole_p", ],
            [ "linear", ],
            ]
    curve_tag_list = [
            ("all", "tt",),
            (1.0, "tt",),
            (1.5, "tt",),
            (2.0, "tt",),
            (2.5, "tt",),
            (10.0, "tt",),
            ]
    t_s_fm_list = [ 0.5, 1.0, 1.5, 2.0, 2.5, ]
    #
    # for pion_type in pion_type_list:
    #     for ff_tag_list in ff_tag_list_list:
    #         for curve_tag in curve_tag_list:
    #             for t_s_fm in t_s_fm_list:
    #                 fit_r_pi_fm(job_tag, jk_list_four_point_em, pion_type, ff_tag_list, curve_tag, t_s_fm)
    pion_type = "pion-diff"
    ff_tag_list = [ "pole", ]
    curve_tag = ("all", "tt",)
    t_s_fm_list = list(np.arange(0.0, 2.6, 0.1))
    table = []
    for t_s_fm in t_s_fm_list:
        jks = fit_r_pi_fm(job_tag, jk_list_four_point_em, pion_type, ff_tag_list, curve_tag, t_s_fm)
        v = [ t_s_fm, q.jk_avg(jks), q.jk_err(jks), ]
        table.append(v)
    content = "\n".join([ " ".join(map(str, v)) for v in table ] + [ "", ])
    q.qtouch_info(get_save_path(f"curve-fits/{job_tag}/curve.txt"), content)
    #
    q.clean_cache()
    q.timer_display()

# PDG: r_pi_fm = 0.659(4)

qg.begin_with_gpt()

q.check_time_limit()

job_tag_list = [
        # "24D",
        # "32D",
        # "24DH",
        # "32Dfine",
        # "48I",
        "64I",
        ]

for job_tag in job_tag_list:
    run_job(job_tag)

qg.end_with_gpt()
