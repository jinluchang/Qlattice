#!/usr/bin/env python3

import os

os.environ['JAX_ENABLE_X64'] = 'True'
os.environ['JAX_PLATFORMS'] = 'cpu'

json_results = []
check_eps = 5e-5

from qlat_utils.q_hlt_reconstruction import *

import numpy as np
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4]]

q.begin_with_mpi(size_node_list)

def mk_f_delta_target(center, sigma):
    def f_delta_target(e):
        return np.exp(-((e - center) / sigma)**2)
    return f_delta_target

f_delta_target = mk_f_delta_target(0.20, 0.15)

t_size = 64
t_arr = np.arange(t_size)
cov = np.eye(t_size) * np.exp(-0.2 * t_arr) * 1e-6

params = mk_hlt_params()
params["f_delta_target"] = f_delta_target
params["t_arr"] = t_arr
params["cov"] = cov
params["alpha"] = -0.01
params["lambda"] = 0.5

g_t_arr_o = mk_g_t_arr(params)

q.json_results_append(f"g_t_arr_o", q.get_data_sig(g_t_arr_o, q.RngState()), check_eps)

aa = aa_from_g(g_t_arr_o, params).item()
aa0 = aa_from_g(np.zeros_like(g_t_arr_o), params).item()
q.displayln_info(f"A={aa} ; A_0={aa0} ; d={np.sqrt(aa / aa0)}")

q.json_results_append(f"aa", aa, check_eps)
q.json_results_append(f"aa0", aa0, check_eps)
q.json_results_append(f"aa / aa0", aa / aa0, check_eps)

e_arr = np.linspace(0.0, 10.0, 1000)

delta_target = f_delta_target(e_arr)

q.json_results_append(f"delta_target", q.get_data_sig(delta_target, q.RngState()), check_eps)

delta = delta_from_g(g_t_arr_o, t_arr, e_arr)
delta = np.array(delta)

q.json_results_append(f"delta", q.get_data_sig(delta, q.RngState()), check_eps)

params = mk_hlt_params()
params["f_delta_target"] = f_delta_target
params["t_arr"] = t_arr
params["cov"] = cov
params["alpha"] = -0.01
params["lambda"] = 0.5
params["e_arr"] = e_arr

g_t_arr_o_via_sum = mk_g_t_arr_via_sum(params)

aa = aa_from_g_via_sum(g_t_arr_o_via_sum, params).item()
aa0 = aa_from_g_via_sum(np.zeros_like(g_t_arr_o_via_sum), params).item()
q.displayln_info(f"A={aa} ; A_0={aa0} ; d={np.sqrt(aa / aa0)} via sum")

q.json_results_append(f"aa ; via sum", aa, check_eps)
q.json_results_append(f"aa0 ; via sum", aa0, check_eps)
q.json_results_append(f"aa / aa0 ; via sum", aa / aa0, check_eps)

aa = aa_from_g(g_t_arr_o_via_sum, params).item()
aa0 = aa_from_g(np.zeros_like(g_t_arr_o_via_sum), params).item()
q.displayln_info(f"A={aa} ; A_0={aa0} ; d={np.sqrt(aa / aa0)} via sum")

q.json_results_append(f"aa ; g via sum", aa, check_eps)
q.json_results_append(f"aa0 ; g via sum", aa0, check_eps)
q.json_results_append(f"aa / aa0 ; g via sum", aa / aa0, check_eps)

delta_via_sum = delta_from_g(g_t_arr_o_via_sum, t_arr, e_arr)
delta_via_sum = np.array(delta_via_sum)

q.json_results_append(f"delta via sum", q.get_data_sig(delta_via_sum, q.RngState()), check_eps)

q.check_log_json(__file__)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
