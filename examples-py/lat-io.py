#!/usr/bin/env python3

import qlat as q
import numpy as np
import os

q.begin()

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

q.displayln_info(f"CHECK: v1 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

ld.save("results/v1.lat")

ld = load_lat_data("results/v1.lat")

q.displayln_info(f"CHECK: v2 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

rs.copy().g_rand_fill(np.asarray(ld))

q.displayln_info(f"CHECK: v3 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

q.displayln_info(str(ld))

ld.save("results/v3.lat")
ld = load_lat_data("results/v3.lat")

q.displayln_info(f"CHECK: v4 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

rs.copy().g_rand_fill(np.asarray(ld).ravel())

q.displayln_info(f"CHECK: v5 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

ld = q.mk_lat_data(info_list, is_complex = False)

q.displayln_info(f"CHECK: v6 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

ld.save("results/v6.lat")
ld = load_lat_data("results/v6.lat")

q.displayln_info(f"CHECK: v7 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

rs.copy().g_rand_fill(np.asarray(ld))

q.displayln_info(f"CHECK: v8 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

ld.save("results/v8.lat")
ld = load_lat_data("results/v8.lat")

q.displayln_info(f"CHECK: v9 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

info_list = []

ld = q.mk_lat_data(info_list, is_complex = False)

q.displayln_info(f"CHECK: v10 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

ld.save("results/v10.lat")
ld = load_lat_data("results/v10.lat")

q.displayln_info(f"CHECK: v11 {ld.info()} {ld.is_complex()} {q.get_double_sig(ld, rs.copy()):.13E}")

q.check_all_files_crc32_info("results")

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end()
