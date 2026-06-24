#!/usr/bin/env python3

import qlat as q
import numpy as np

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")

rs = q.RngState("test")

info_list = [
    [
        "dim1",
        4,
    ],
    [
        "dim2 name",
        3,
    ],
    [
        "dim3",
        2,
        [
            "u",
            "d",
        ],
    ],
    [
        "dim4",
        5,
    ],
]

ld = q.mk_lat_data(info_list)
rs.copy().g_rand_fill(np.asarray(ld).ravel())
q.json_results_append("lat-io-test2: v0 qnorm", ld.qnorm())
q.json_results_append(f"lat-io-test2: v0 ndim={ld.ndim()}")
q.json_results_append(f"lat-io-test2: v0 dim_sizes={ld.dim_sizes()}")
q.json_results_append(f"lat-io-test2: v0 info={ld.info()} is_complex={ld.is_complex()}")
q.json_results_append("lat-io-test2: v0 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
q.json_results_append(
    "lat-io-test2: v0 data_sig rs1", q.get_data_sig(ld, rs.split("1")), 1e-14
)
ld.set_zero()
q.json_results_append(f"lat-io-test2: v1 info={ld.info()} is_complex={ld.is_complex()}")
q.json_results_append("lat-io-test2: v1 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
ld.save("results/v1.lat")
ld = q.load_lat_data("results/v1.lat")
q.json_results_append(f"lat-io-test2: v2 info={ld.info()} is_complex={ld.is_complex()}")
q.json_results_append("lat-io-test2: v2 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
rs.copy().g_rand_fill(np.asarray(ld))
q.json_results_append("lat-io-test2: v3 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
q.json_results_append(f"lat-io-test2: v3 info={ld.info()} is_complex={ld.is_complex()}")
ld.save("results/v3.lat")
ld = q.load_lat_data("results/v3.lat")
q.json_results_append(f"lat-io-test2: v4 info={ld.info()} is_complex={ld.is_complex()}")
q.json_results_append("lat-io-test2: v4 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
rs.copy().g_rand_fill(np.asarray(ld).ravel())
arr = ld.to_numpy()
ld.from_numpy(arr)
q.json_results_append("lat-io-test2: v5 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
q.json_results_append(f"lat-io-test2: v5 info={ld.info()} is_complex={ld.is_complex()}")

ld = q.mk_lat_data(info_list, is_complex=False)
q.json_results_append(f"lat-io-test2: v6 info={ld.info()} is_complex={ld.is_complex()}")
q.json_results_append("lat-io-test2: v6 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
ld.save("results/v6.lat")
ld = q.load_lat_data("results/v6.lat")
q.json_results_append(f"lat-io-test2: v7 info={ld.info()} is_complex={ld.is_complex()}")
q.json_results_append("lat-io-test2: v7 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
rs.copy().g_rand_fill(np.asarray(ld))
q.json_results_append(f"lat-io-test2: v8 info={ld.info()} is_complex={ld.is_complex()}")
q.json_results_append("lat-io-test2: v8 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
ld.save("results/v8.lat")
ld = q.load_lat_data("results/v8.lat")
arr = ld.to_numpy()
ld.from_numpy(arr, is_complex=False)
q.json_results_append("lat-io-test2: v9 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14)
q.json_results_append(f"lat-io-test2: v9 info={ld.info()} is_complex={ld.is_complex()}")

info_list = [
    [
        "dim",
        3,
        [
            "a",
            "b",
            "c",
        ],
    ],
]

ld = q.mk_lat_data(info_list)
q.json_results_append(
    f"lat-io-test2: v10 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v10 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)
rs.copy().g_rand_fill(np.asarray(ld))
q.json_results_append(
    f"lat-io-test2: v11 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v11 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)
ld.save("results/v11.lat")
ld = q.load_lat_data("results/v11.lat")
arr = ld.to_numpy()
ld.from_numpy(arr)
q.json_results_append(
    f"lat-io-test2: v12 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v12 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)

ld = q.mk_lat_data(info_list, is_complex=False)
q.json_results_append(
    f"lat-io-test2: v13 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v13 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)
ld.save("results/v13.lat")
ld = q.load_lat_data("results/v13.lat")
arr = ld.to_numpy()
ld.from_numpy(arr, is_complex=False)
q.json_results_append(
    f"lat-io-test2: v14 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v14 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)

info_list = []

ld = q.mk_lat_data(info_list)
q.json_results_append(
    f"lat-io-test2: v15 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v15 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)
rs.copy().g_rand_fill(np.asarray(ld))
q.json_results_append(
    f"lat-io-test2: v16 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v16 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)
ld.save("results/v16.lat")
ld = q.load_lat_data("results/v16.lat")
arr = ld.to_numpy()
ld.from_numpy(arr)
q.json_results_append(
    f"lat-io-test2: v17 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v17 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)

ld = q.mk_lat_data(info_list, is_complex=False)
q.json_results_append(
    f"lat-io-test2: v18 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v18 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)
rs.copy().g_rand_fill(np.asarray(ld))
q.json_results_append(
    f"lat-io-test2: v19 info={ld.info()} is_complex={ld.is_complex()}"
)
q.json_results_append(
    "lat-io-test2: v19 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)
ld.save("results/v19.lat")
ld = q.load_lat_data("results/v19.lat")
arr = ld.to_numpy()
ld.from_numpy(arr, is_complex=False)
q.json_results_append(
    "lat-io-test2: v20 data_sig", q.get_data_sig(ld, rs.copy()), 1e-14
)
q.json_results_append(
    f"lat-io-test2: v20 info={ld.info()} is_complex={ld.is_complex()}"
)

q.check_all_files_crc32_info("results")

q.check_log_json(__file__, check_eps=1e-14)
q.timer_display()

q.end_with_mpi()

q.displayln_info("CHECK: finished successfully.")
