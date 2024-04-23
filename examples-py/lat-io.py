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

ld = q.LatData()

ld.load("results/ld.lat")

assert q.qnorm(np.asarray(ld) - arr) == 0

ld1 = q.load_lat_data("results/ld.lat")

assert q.qnorm(np.asarray(ld1) - arr) == 0

ld1 = q.mk_lat_data_real_f(info_list)

q.displayln_info(ld1.info())

ld1 @= ld

arr1 = ld1.to_numpy()

ld1.save("results/ld1.latf")

ld1 = q.LatDataRealF()

ld1.load("results/ld1.latf")

assert q.qnorm(np.asarray(ld1) - arr1) == 0

ld2 = q.mk_lat_data_long(info_list)
arr2 = ld2[:]
arr2.ravel()[:] = np.arange(len(arr2.ravel()))

ld2.save("results/ld2.latl")

ld2 = q.LatDataLong()
ld2.load("results/ld2.latl")

assert q.qnorm(ld2[:] - arr2) == 0

ld3 = q.mk_lat_data_int(info_list)
arr3 = ld3[:]
arr3.ravel()[:] = np.arange(len(arr3.ravel()))

ld3.save("results/ld3.latl")

ld3 = q.LatDataInt()
ld3.load("results/ld3.latl")

assert q.qnorm(ld3[:] - arr3) == 0

q.check_all_files_crc32_info("results")

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
