#!/usr/bin/env python3

import qlat as q
import numpy as np

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin_with_mpi(size_node_list)

q.json_results_append(q.show_val(1))
q.json_results_append(q.show_val(1132))
q.json_results_append(q.show_val(11.32))
q.json_results_append(q.show_val(-11.32))
q.json_results_append(q.show_val(-11.32, num_float_digit=4))
q.json_results_append(q.show_val(-11.32, num_float_digit=1))
q.json_results_append(q.show_val(-11.32, num_float_digit=True))
q.json_results_append(q.show_val(-11.32, num_exp_digit=4))
q.json_results_append(q.show_val(-11.32, num_exp_digit=1))
q.json_results_append(q.show_val(-11.32, num_exp_digit=True))
q.json_results_append(q.show_val(-11.32, exponent=3))
q.json_results_append(q.show_val(-11.32, exponent=3, is_latex=False))

q.json_results_append(q.show_val_err(1))
q.json_results_append(q.show_val_err(1132))
q.json_results_append(q.show_val_err(11.32))
q.json_results_append(q.show_val_err((11.32, 12)))
q.json_results_append(q.show_val_err((-11.32, 120)))
q.json_results_append(q.show_val_err((1.12e16, 12), num_float_digit=1))
q.json_results_append(q.show_val_err((1.12e16, 12e6), num_exp_digit=True))
q.json_results_append(q.show_val_err((1.12e16, 12e6)))
q.json_results_append(q.show_val_err((1.12e16, 12e7), exponent=10))
q.json_results_append(q.show_val_err((1.12e16, 12e7), exponent=10, is_latex=False))

q.timer_display()
q.check_log_json(__file__, check_eps=1e-5)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
