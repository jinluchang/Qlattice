#!/usr/bin/env python3

import qlat as q
import numpy as np
import os

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")

q.benchmark_matrix_functions(16 * 1024)

q.check_all_files_crc32_info("results")

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
