#!/usr/bin/env python3

import qlat as q
import numpy as np

q.begin_with_mpi()

data_x_arr = np.linspace(0.2, 3.4, 20)
data_arr = np.array([ np.sin(data_x_arr), np.cos(data_x_arr), ])

x_arr = np.linspace(0.0, 4.0, 10)
i_arr = q.interp_i_arr(data_x_arr, x_arr)
approx_x_arr = q.interp(data_x_arr, i_arr)

q.json_results_append(f"i_arr", i_arr)
q.json_results_append(f"approx_x_arr", approx_x_arr)

interp_data_arr = q.interp(data_arr, i_arr)

q.json_results_append(f"interp_data_arr", interp_data_arr.ravel())

interp2_data_arr = q.interp_x(data_arr, data_x_arr, x_arr)
assert np.all(interp2_data_arr == interp_data_arr)

data_arr = data_arr.T
interp_data_arr = q.interp(data_arr, i_arr, axis=0)
q.json_results_append(f"interp_data_arr", interp_data_arr.ravel())

interp2_data_arr = q.interp_x(data_arr, data_x_arr, x_arr, axis=0)
assert np.all(interp2_data_arr == interp_data_arr)

q.check_log_json(__file__, check_eps=1e-10)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
