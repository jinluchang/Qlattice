#!/usr/bin/env python3

import os

os.environ['JAX_ENABLE_X64'] = 'True'
os.environ['JAX_PLATFORMS'] = 'cpu'

json_results = []
check_eps = 5e-5

def json_results_append(*args):
    q.displayln_info(r"//------------------------------------------------------------\\")
    q.displayln_info(-1, *args)
    q.displayln_info(r"\\------------------------------------------------------------//")
    json_results.append(args)

import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4]]

q.begin_with_mpi(size_node_list)

mp_pool_n_proc = q.get_q_num_mp_processes()
if mp_pool_n_proc > 0:
    import multiprocessing
    mp_pool = multiprocessing.get_context('fork').Pool(mp_pool_n_proc, initializer=q.q_fit_corr_2.mp_initializer)
else:
    mp_pool = None

q.displayln_info(f"mp_pool_n_proc={mp_pool_n_proc}")

n_eigs = 4

n_ops = 3

t_size = 10

param_arr_setup, jk_corr_data, corr_data_sigma, t_arr = q.q_fit_corr_2.mk_data_set(
    n_jk=8, n_ops=n_ops, n_eigs=n_eigs, t_size=t_size,
    sigma=1e-3,
    rng=q.RngState("seed-data-set-4"))

param_arr_setup = q.q_fit_corr_2.sort_param_arr_free_eig(param_arr_setup, n_ops, free_eig_idx_arr=np.arange(n_eigs))

json_results_append(f"n_eigs = {n_eigs}")

for i in range(n_eigs):
    json_results_append(f"param_arr_setup[{i}]", param_arr_setup[i], 1e-7)

sig1 = q.get_data_sig(param_arr_setup, q.RngState())

json_results_append(f"param_arr_setup sig1", sig1, 1e-7)

sig2 = q.get_data_sig(jk_corr_data, q.RngState())

json_results_append(f"jk_corr_data sig2", sig2, 1e-7)

sig3 = q.get_data_sig(corr_data_sigma, q.RngState())

json_results_append(f"corr_data_sigma sig3", sig3, 1e-7)

e_arr = param_arr_setup[:n_eigs]

res = q.q_fit_corr_2.fit_eig_coef(
        jk_corr_data[:, :, :, 1:10],
        t_arr=t_arr[1:10],
        e_arr=e_arr,
        t_size=t_size,
        free_eig_idx_arr=[ 1, ],
        n_step_mini_avg=5,
        n_step_mini_jk=5,
        )

sig4 = q.get_data_sig(res['jk_chisq'], q.RngState())

json_results_append(f"jk_chisq sig4", sig4, 1e-4)

sig5 = q.get_data_sig(abs(res['jk_param_arr']), q.RngState())

json_results_append(f"jk_param_arr sig5", sig5, 1e-4)

sig6 = q.get_data_sig(abs(res['jk_chisq_grad']), q.RngState())

json_results_append(f"jk_chisq_grad sig6", sig6, 1e-4)

sig7 = q.get_data_sig(abs(res['jk_param_arr_for_scaled_corr']), q.RngState())

json_results_append(f"jk_param_arr_for_scaled_corr sig7", sig7, 1e-4)

param_compare = np.stack((param_arr_setup,) + q.g_jk_avg_err(res['jk_param_arr'])).T

q.displayln_info(f"{param_compare}")

res = q.q_fit_corr_2.fit_eig_coef(
        jk_corr_data[:, :, :, 1:10],
        t_arr=t_arr[1:10],
        e_arr=e_arr,
        t_size=t_size,
        free_eig_idx_arr=[ 1, ],
        n_step_mini_avg=5,
        n_step_mini_jk=5,
        mp_pool=mp_pool,
        )

sig4 = q.get_data_sig(res['jk_chisq'], q.RngState())

json_results_append(f"jk_chisq sig4", sig4, 1e-4)

sig5 = q.get_data_sig(abs(res['jk_param_arr']), q.RngState())

json_results_append(f"jk_param_arr sig5", sig5, 1e-4)

sig6 = q.get_data_sig(abs(res['jk_chisq_grad']), q.RngState())

json_results_append(f"jk_chisq_grad sig6", sig6, 1e-4)

sig7 = q.get_data_sig(abs(res['jk_param_arr_for_scaled_corr']), q.RngState())

json_results_append(f"jk_param_arr_for_scaled_corr sig7", sig7, 1e-4)

res = q.q_fit_corr_2.fit_eig_coef(
        jk_corr_data[:, :, :, 1:10],
        t_arr=t_arr[1:10],
        e_arr=e_arr,
        t_size=t_size,
        free_eig_idx_arr=[ 1, ],
        eig_maximum_arr=[ 0.5, ],
        n_step_mini_avg=5,
        n_step_mini_jk=5,
        mp_pool=mp_pool,
        )

sig8 = q.get_data_sig(res['jk_chisq'], q.RngState())

json_results_append(f"jk_chisq sig8", sig8, 1e-2)

sig9 = q.get_data_sig(abs(res['jk_param_arr']), q.RngState())

json_results_append(f"jk_param_arr sig9", sig9, 1e-4)

q.check_log_json(__file__, json_results)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
