#!/usr/bin/env python3

import qlat as q

q.begin_with_mpi()

q.benchmark_matrix_functions(16 * 1024)

q.json_results_append("benchmark-matrix-functions")
q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
