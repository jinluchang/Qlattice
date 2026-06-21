#!/usr/bin/env python3

import qlat as q

q.begin_with_mpi()

q.show_machine()

q.json_results_append('check-id-node')

q.check_log_json(__file__)

q.end_with_mpi()

q.displayln_info("CHECK: finished successfully.")
