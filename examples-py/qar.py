#!/usr/bin/env python3

import qlat as q
import os
import numpy as np

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")

content = b"hello world!\n"

if q.get_id_node() == 0:
    qfile = q.open_qfile("results/test-qfile.txt", "w")
    qfile.write(content)
    qfile.close()

assert content == q.qcat_bytes_sync_node("results/test-qfile.txt")

if q.get_id_node() == 0:
    qfile0 = q.open_qfile("results/test-qfile.txt", "r")
    qfile = q.open_qfile("results/test-qfile-2.txt", "w")
    qfile.write(qfile0)
    qfile0.close()
    qfile.close()

assert content == q.qcat_bytes_sync_node("results/test-qfile-2.txt")

qfile = q.open_qfile_str("/ string /", "w")
qfile.write(content)
assert content == qfile.content_bytes()
qfile.close()

def test_ld_str_io(path):
    ld_load = q.load_lat_data(path)
    ld_bytes = ld_load.save_str()
    assert q.qcat_bytes_sync_node(path) == ld_bytes
    ld_str = q.LatData()
    ld_str.load_str(ld_bytes)
    assert ld_str.save_str() == ld_bytes

ld = q.LatData()
ld.save(f"results/data/ld.lat")

test_ld_str_io(f"results/data/ld.lat")

for i in range(20):
    ld = q.LatData()
    rs = q.RngState(f"seed {i}")
    arr = rs.g_rand_arr((2, 5, 10, 10,))
    ld.from_numpy(arr)
    ld.save(f"results/data/ld-1000/ld-{i}-1000.lat")
    test_ld_str_io(f"results/data/ld-1000/ld-{i}-1000.lat")

ld = q.LatData()
ld.from_numpy(np.arange(10000.0).astype(complex).reshape(2, 5, 10, 100))
ld.save(f"results/data/ld-10000.lat")
test_ld_str_io(f"results/data/ld-10000.lat")

def test_ldi_str_io(path):
    ld_load = q.load_lat_data_int(path)
    ld_bytes = ld_load.save_str()
    assert q.qcat_bytes_sync_node(path) == ld_bytes
    ld_str = q.LatDataInt()
    ld_str.load_str(ld_bytes)
    assert ld_str.save_str() == ld_bytes

ld = q.LatDataInt()
ld.from_numpy(np.arange(10000).reshape(2, 5, 10, 100))
ld.save(f"results/data/ld-10000.lati")
test_ldi_str_io(f"results/data/ld-10000.lati")

if q.get_id_node() == 0:
    qar = q.open_qar("results/test-qar.qar", "w")
    qar.write("f1", "f1-info", content)
    qar.write("f2", "f2-info", content)
    l0 = qar.list()
    q.displayln_info(f"CHECK: l0 {l0}")
    qar.close()
    l1 = q.list_qar("results/test-qar.qar")
    q.displayln_info(f"CHECK: l1 {l1}")
    assert l0 == l1
    qar = q.open_qar("results/test-qar.qar", "w")
    qar.write("f1", "f1-info", content)
    qar.close()
    l2 = q.list_qar("results/test-qar.qar")
    q.displayln_info(f"CHECK: l2 {l2}")
    qar = q.open_qar("results/test-qar.qar", "a")
    qar.write("f2", "f2-info", content)
    l3 = qar.list()
    q.displayln_info(f"CHECK: l3 {l3}")
    qar.close()
    l4 = q.list_qar("results/test-qar.qar")
    q.displayln_info(f"CHECK: l4 {l4}")
    qfile = q.open_qfile("results/test-qar.qar", "a")
    qfile.write("hello")
    qfile.close()
    qar = q.open_qar("results/test-qar.qar", "a")
    qar.write("f3", "f3-info", content)
    l5 = qar.list()
    q.displayln_info(f"CHECK: l5 {l5}")
    qar.close()
    l6 = q.list_qar("results/test-qar.qar")
    q.displayln_info(f"CHECK: l6 {l6}")
    assert l5 == l6

qar = q.open_qar_info("results/test-qar.qar", "r")
l = qar.list()
l = q.bcast_py(l)
q.displayln_info(f"CHECK: open_qar_info l='{l}'")
for fn in l:
    data = qar.read_data(fn)
    q.displayln_info(f"CHECK: open_qar_info fn='{fn}' data='{data!r}'")
qar.close()

q.qar_create_info(f"results/data.qar", f"results/data")

q.qar_extract_info(f"results/data.qar", f"results/data2")

q.qar_create_info(f"results/data2.qar", f"results/data2", is_remove_folder_after=True)

q.qar_extract_info(f"results/data2.qar", f"results/data2", is_remove_qar_after=True)

qar_multi_vol_max_size = q.get_qar_multi_vol_max_size()
q.displayln_info(f"CHECK: qar_multi_vol_max_size={qar_multi_vol_max_size}")

q.set_qar_multi_vol_max_size(16 * 1024)
qar_multi_vol_max_size = q.get_qar_multi_vol_max_size()
q.displayln_info(f"CHECK: qar_multi_vol_max_size={qar_multi_vol_max_size}")

q.qar_create_info(f"results/data2/ld-1000.qar", f"results/data2/ld-1000", is_remove_folder_after=True)

q.qar_create_info(f"results/data2.qar", f"results/data2")

q.qar_extract_info(f"results/data2.qar", f"results/data3")

q.qar_create_info(f"results/data3.qar", f"results/data3", is_remove_folder_after=True)

q.qar_extract_info(f"results/data3.qar", f"results/data3", is_remove_qar_after=True)

q.qar_create_info(f"results/data4.qar", f"results/data3")

ld = q.LatData()
ld.load(f"results/data4/ld-10000.lat")
assert q.qnorm(q.load_lat_data(f"results/data4/ld-10000.lat") - q.load_lat_data(f"results/data/ld-10000.lat")) == 0

l1 = q.list_qar("results/data4.qar")

q.sync_node()
q.displayln_info("CHECK: ", l1)
q.sync_node()

l2 = q.list_qar("results/data4/ld-1000.qar")

q.sync_node()
q.displayln_info("CHECK: ", l2)
q.sync_node()

num_clean_up_qfiles = q.clean_up_qfile_map()
q.displayln_info(f"CHECK: clean_up_qfile_map num_clean_up_qfiles={num_clean_up_qfiles}");

sq_list = sorted(q.show_all_qfile())
q.sync_node()
q.displayln_info(f"CHECK: q.show_all_qfile()")
for idx, s in enumerate(sq_list):
    q.displayln_info(f"CHECK: {idx} {s}")
q.sync_node()

for fn in [ f"ld-10000.lat", f"ld-1000/ld-1-1000.lat", ]:
    q.qcopy_file_info(f"results/data4/{fn}", f"results/data5/{fn}")
    if 0 == q.get_id_node():
        assert q.qcat_bytes(f"results/data4/{fn}") == q.qcat_bytes(f"results/data5/{fn}")

q.check_all_files_crc32_info("results")

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
