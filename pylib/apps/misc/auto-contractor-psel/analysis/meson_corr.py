#!/usr/bin/env python3

import qlat as q
import rbc_ukqcd_params as rup
import os
import numpy as np

save_path_default = "results"

load_path_list = [
        "results",
        "../results",
        "/home/frank/qcdqedta/luchang/auto-contractor-psel/v5/results",
        ]

def get_save_path(fn):
    return os.path.join(save_path_default, fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = load_path_list
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

def check_traj(job_tag, traj):
    # return True is all files are available
    fn_list = [
            f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr/wsnk_wsrc.lat",
            f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr/psnk_wsrc.lat",
            f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr_with_env/wsnk_wsrc.lat",
            f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr_with_env/psnk_wsrc.lat",
            ]
    for fn in fn_list:
        if get_load_path(fn) is None:
            # q.displayln_info(f"check_traj: {fn} does not exist.")
            return False
    return True

@q.timer_verbose
def get_trajs(job_tag):
    all_trajs = range(200, 3000, 5)
    trajs = []
    for traj in all_trajs:
        if check_traj(job_tag, traj):
            trajs.append(traj)
    return trajs

@q.timer
def get_data(job_tag, traj):
    ld_ww = q.load_lat_data(get_load_path(f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr_with_env/wsnk_wsrc.lat"))
    ld_pw = q.load_lat_data(get_load_path(f"auto-contractor-psel/{job_tag}/traj={traj}/meson_corr_with_env/psnk_wsrc.lat"))
    return [ ld_ww, ld_pw, ]

def get_all_data(job_tag, trajs):
    # lds_all = [ lds_ww, lds_pw, ]
    # lds_ww = [ ld_ww, ... ]
    # lds_pw = [ ld_pw, ... ]
    # ld_ww[i_k][tsep] = val
    lds_all = [ [], [], ]
    for traj in trajs:
        for lds, ld in zip(lds_all, get_data(job_tag, traj)):
            lds.append(q.Data(ld[(0,)][:,:,0]))
    lds_all = [ q.jackknife(lds) for lds in lds_all ]
    return lds_all

def get_tsep_env(job_tag):
    return 6

def subtract_discon(job_tag, ld_ww, ld_pw):
    tsep_env = get_tsep_env(job_tag)
    t_size = ld_ww.shape[1]
    n_inner = 4
    i_inner_pipi = 0
    i_inner_apipi = 1
    i_inner_kk = 2
    i_inner_akk = 3
    i_outer_same = 1
    i_outer_other = 2
    i_outer_oppo = 3
    for tsep in range(t_size):
        t_env = (tsep + 2 * tsep_env) % t_size
        t_env1 = (tsep + tsep_env) % t_size
        ld_ww[i_outer_same * n_inner + i_inner_pipi][tsep] -= (
                ld_ww[i_inner_pipi][tsep] * ld_ww[i_inner_pipi][t_env]
                + ld_ww[i_inner_pipi][t_env1] * ld_ww[i_inner_pipi][t_env1])
        ld_ww[i_outer_other * n_inner + i_inner_pipi][tsep] -= (
                ld_ww[i_inner_pipi][tsep] * ld_ww[i_inner_pipi][t_env])
        ld_ww[i_outer_oppo * n_inner + i_inner_pipi][tsep] -= (
                ld_ww[i_inner_pipi][tsep] * ld_ww[i_inner_pipi][t_env]
                + ld_ww[i_inner_pipi][tsep_env] * ld_ww[i_inner_pipi][tsep_env])
        ld_ww[i_outer_same * n_inner + i_inner_apipi][tsep] -= (
                ld_ww[i_inner_apipi][tsep] * ld_ww[i_inner_pipi][t_env]
                + ld_ww[i_inner_apipi][t_env1] * ld_ww[i_inner_pipi][t_env1])
        ld_ww[i_outer_other * n_inner + i_inner_apipi][tsep] -= (
                ld_ww[i_inner_apipi][tsep] * ld_ww[i_inner_pipi][t_env])
        ld_ww[i_outer_oppo * n_inner + i_inner_apipi][tsep] -= (
                ld_ww[i_inner_apipi][tsep] * ld_ww[i_inner_pipi][t_env]
                - ld_ww[i_inner_apipi][tsep_env] * ld_ww[i_inner_pipi][tsep_env])
        ld_ww[i_outer_same * n_inner + i_inner_pipi][tsep] -= (
                ld_ww[i_inner_pipi][tsep] * ld_ww[i_inner_pipi][t_env]
                + ld_ww[i_inner_pipi][t_env1] * ld_ww[i_inner_pipi][t_env1])
        ld_pw[i_outer_other * n_inner + i_inner_pipi][tsep] -= (
                ld_pw[i_inner_pipi][tsep] * ld_ww[i_inner_pipi][t_env])
        ld_pw[i_outer_oppo * n_inner + i_inner_pipi][tsep] -= (
                ld_pw[i_inner_pipi][tsep] * ld_ww[i_inner_pipi][t_env]
                + ld_pw[i_inner_pipi][tsep_env] * ld_ww[i_inner_pipi][tsep_env])
        ld_pw[i_outer_same * n_inner + i_inner_apipi][tsep] -= (
                ld_pw[i_inner_apipi][tsep] * ld_ww[i_inner_pipi][t_env]
                + ld_pw[i_inner_apipi][t_env1] * ld_ww[i_inner_pipi][t_env1])
        ld_pw[i_outer_other * n_inner + i_inner_apipi][tsep] -= (
                ld_pw[i_inner_apipi][tsep] * ld_ww[i_inner_pipi][t_env])
        ld_pw[i_outer_oppo * n_inner + i_inner_apipi][tsep] -= (
                ld_pw[i_inner_apipi][tsep] * ld_ww[i_inner_pipi][t_env]
                - ld_pw[i_inner_apipi][tsep_env] * ld_ww[i_inner_pipi][tsep_env])
        for i_inner in [ i_inner_kk, i_inner_akk, ]:
            for i_outer in [ 1, 2, 3, ]:
                ld_ww[i_outer * n_inner + i_inner][tsep] -= (
                        ld_ww[i_inner][tsep] * ld_ww[i_inner_pipi][t_env])
                ld_pw[i_outer * n_inner + i_inner][tsep] -= (
                        ld_pw[i_inner][tsep] * ld_ww[i_inner_pipi][t_env])

def divide_norm(job_tag, ld_ww, ld_pw):
    tsep_env = get_tsep_env(job_tag)
    t_size = ld_ww.shape[1]
    n_inner = 4
    i_inner_pipi = 0
    for i_outer in [ 1, 2, 3, ]:
        for i_inner in [ 0, 1, 2, 3, ]:
            for tsep in range(t_size):
                t_env = (tsep + 2 * tsep_env) % t_size
                ld_ww[i_outer * n_inner + i_inner][tsep] /= ld_ww[i_inner_pipi][t_env]
                ld_pw[i_outer * n_inner + i_inner][tsep] /= ld_ww[i_inner_pipi][t_env]
                ld_ww[i_outer * n_inner + i_inner][tsep] /= ld_ww[i_inner][tsep]
                ld_pw[i_outer * n_inner + i_inner][tsep] /= ld_pw[i_inner][tsep]

def analysis(job_tag):
    trajs = get_trajs(job_tag)
    q.displayln_info(f"analysis: trajs={trajs}")
    if not trajs:
        return
    lds_all = get_all_data(job_tag, trajs)
    [ lds_ww, lds_pw, ] = lds_all
    for ld_ww, ld_pw in zip(lds_ww, lds_pw):
        subtract_discon(job_tag, ld_ww.val, ld_pw.val)
        divide_norm(job_tag, ld_ww.val, ld_pw.val)
        pass
    data = []
    i_inner = 1
    for ld_ww, ld_pw in zip(lds_ww, lds_pw):
        ld = sum([ ld_pw.val[i * 4 + i_inner] for i in [ 1, 2, 3, ] ])
        data.append(ld)
    q.displayln_info(q.jk_avg(lds_pw).val[i_inner,:20])
    q.displayln_info(q.jk_avg(data)[:20])
    q.displayln_info(q.jk_err(data)[:20])

job_tag = "24D"

analysis(job_tag)
