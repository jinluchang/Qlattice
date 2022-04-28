#!/usr/bin/env python3

import qlat as q
import os

q.begin()

q.qremove_all_info("results")
q.qmkdir_info("results")

ld = q.LatData()

ld.save(f"results/ld.lat")

q.check_all_files_crc32_info("results")

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end()
