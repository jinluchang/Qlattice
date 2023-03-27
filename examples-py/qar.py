#!/usr/bin/env python3

import qlat as q
import os
import numpy as np

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")

ld = q.LatData()
ld.save(f"results/data/ld.lat")

for i in range(20):
    ld = q.LatData()
    ld.from_numpy(np.arange(1000.0).astype(complex).reshape(2, 5, 10, 10))
    ld.save(f"results/data/ld-1000/ld-{i}-1000.lat")

ld = q.LatData()
ld.from_numpy(np.arange(10000.0).astype(complex).reshape(2, 5, 10, 100))
ld.save(f"results/data/ld-10000.lat")

q.qar_create_info(f"results/data.qar", f"results/data")

q.qar_extract_info(f"results/data.qar", f"results/data2")

q.qar_create_info(f"results/data2.qar", f"results/data2", is_remove_folder_after = True)

q.qar_extract_info(f"results/data2.qar", f"results/data2", is_remove_qar_after = True)

qar_multi_vol_max_size = q.get_qar_multi_vol_max_size(16 * 1024)
q.displayln_info(f"CHECK: qar_multi_vol_max_size={qar_multi_vol_max_size}")

q.qar_create_info(f"results/data2/ld-1000.qar", f"results/data2/ld-1000", is_remove_folder_after = True)

q.qar_create_info(f"results/data2.qar", f"results/data2")

q.qar_extract_info(f"results/data2.qar", f"results/data3")

q.qar_create_info(f"results/data3.qar", f"results/data3", is_remove_folder_after = True)

q.qar_extract_info(f"results/data3.qar", f"results/data3", is_remove_qar_after = True)

q.qar_create_info(f"results/data4.qar", f"results/data3")

ld = q.LatData()
ld.load(f"results/data4/ld-10000.lat")
assert q.qnorm(q.load_lat_data(f"results/data4/ld-10000.lat") - q.load_lat_data(f"results/data/ld-10000.lat")) == 0

q.displayln_info("CHECK: ", q.list_qar("results/data4.qar"))
q.displayln_info("CHECK: ", q.list_qar("results/data4/ld-1000.qar"))

for fn in [ f"ld-10000.lat", f"ld-1000/ld-1-1000.lat", ]:
    q.qcopy_file_info(f"results/data4/{fn}", f"results/data5/{fn}")
    if 0 == q.get_id_node():
        assert q.qcat_bytes(f"results/data4/{fn}") == q.qcat_bytes(f"results/data5/{fn}")

q.check_all_files_crc32_info("results")

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end_with_mpi()
