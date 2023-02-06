#!/usr/bin/env python3

import qlat as q
import os

q.begin()

q.qremove_all_info("results")
q.qmkdir_info("results")

ld = q.LatData()

dim_sizes = [ 2, 4, 3, ]

ld.set_dim_sizes(dim_sizes)
ld.set_dim_name(0, "index0", [ "u", "d" ])
ld.set_dim_name(1, "index1")
ld.set_dim_name(2, "index2")

ld.set_zero()

q.displayln_info(f"CHECK: ld: ndim {ld.ndim()}")
q.displayln_info(f"CHECK: ld: dim_sizes {ld.dim_sizes()}")
for dim in range(ld.ndim()):
    q.displayln_info(f"CHECK: ld: dim_name {dim} {ld.dim_name(dim)} {ld.dim_indices(dim)}")

for i0 in range(dim_sizes[0]):
    for i1 in range(dim_sizes[1]):
        for i2 in range(dim_sizes[2]):
            ld[(i0, i1, i2,)] = 1e6 + i0 * 1e4 + i1 * 1e2 + i2

q.displayln_info("ld:")
q.displayln_info(ld.show())
q.displayln_info(f"CHECK: qnorm = {ld.qnorm()}")

ld[(0, 1,)] = [ i * 3 + i * 1j for i in range(3) ]

q.displayln_info("CHECK: ", ld[(0,)])
q.displayln_info("CHECK: ", ld[(1, 2,)])

ld.save("results/test.lat")

ld = q.LatData()
ld.load("results/test.lat")
q.displayln_info("CHECK: ", ld[(1, 3,)])

ld1 = ld.copy()
ld1 += ld1
ld1 *= 0.5
ld1 -= ld
q.displayln_info(ld1.show())
q.displayln_info(f"CHECK: qnorm = {ld1.qnorm()}")

q.displayln_info(ld.to_list())

ld1.from_list(ld.to_list())

q.displayln_info(ld1.show())
q.displayln_info(f"CHECK: qnorm = {ld1.qnorm()}")

q.check_all_files_crc32_info("results")

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end()
