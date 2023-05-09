#!/usr/bin/env python3

import qlat as q
import numpy as np
import os

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")

rs = q.RngState("test")

info_list = [
        [ "dim1", 4, ],
        [ "dim2 name", 3, ],
        [ "dim3", 2, [ "u", "d", ], ],
        [ "dim4", 5, ],
        ]

ld = q.mk_lat_data(info_list)

q.displayln_info(ld.info())

rs.copy().g_rand_fill(np.asarray(ld).ravel())

arr = ld.to_numpy()

ld.save("results/ld.lat")

ld.set_zero()

ld.load("results/ld.lat")

assert q.qnorm(np.asarray(ld) - arr) == 0

ld1 = q.load_lat_data("results/ld.lat")

assert q.qnorm(np.asarray(ld1) - arr) == 0

q.check_all_files_crc32_info("results")

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end_with_mpi()
