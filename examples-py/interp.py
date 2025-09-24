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

data_arr = np.array([ data_x_arr, data_x_arr**2, data_x_arr**3, ])
threshold = 2.1
i_arr = q.get_threshold_i_arr(data_arr, threshold)
q.json_results_append(f"i_arr", i_arr)
assert np.all(q.get_threshold_i_arr(data_arr.T, threshold, axis=0) == i_arr)

threshold_arr = [ 2.1, 2.2, 2.3, ]
i2_arr = q.get_threshold_i_arr(data_arr, threshold_arr)
q.json_results_append(f"i2_arr", i2_arr)
assert np.all(np.array([ q.get_threshold_idx(data_arr[i], threshold_arr[i]) for i in range(len(threshold_arr)) ], dtype=np.float64) == i2_arr)
assert np.all(q.get_threshold_i_arr(data_arr.T, threshold_arr, axis=0) == i2_arr)

assert np.allclose(np.array([ q.interp(data_arr[i], i2_arr[i]) for i in range(len(i2_arr)) ], dtype=np.float64), threshold_arr, rtol=1e-08, atol=1e-08)

x_arr = q.get_threshold_x_arr(data_arr, data_x_arr, threshold)
q.json_results_append(f"x_arr", x_arr)
assert np.allclose(q.interp(data_x_arr, i_arr), x_arr, rtol=1e-08, atol=1e-08)

x2_arr = q.get_threshold_x_arr(data_arr, data_x_arr, threshold_arr)
q.json_results_append(f"x2_arr", x2_arr)
assert np.allclose(q.interp(data_x_arr, i2_arr), x2_arr, rtol=1e-08, atol=1e-08)

q.check_log_json(__file__, check_eps=1e-10)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
