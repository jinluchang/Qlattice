#!/usr/bin/env python3

import copy
import pickle

import numpy as np
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

q.json_results_append("test lat_data documentation examples")

q.json_results_append("construction")

ld = q.mk_lat_data([["t", 8], ["p", 3]])
assert ld.ndim() == 2
assert ld.dim_sizes() == [8, 3]
q.json_results_append(f"mk_lat_data: ndim={ld.ndim()} dim_sizes={ld.dim_sizes()}")

ld_rf = q.mk_lat_data_real_f([["t", 4], ["x", 2]])
assert ld_rf.ndim() == 2
q.json_results_append(f"mk_lat_data_real_f: ndim={ld_rf.ndim()}")

ld_int = q.mk_lat_data_int([["a", 3]])
assert ld_int.ndim() == 1
assert not ld_int.is_complex()
q.json_results_append(
    f"mk_lat_data_int: ndim={ld_int.ndim()} is_complex={ld_int.is_complex()}"
)

ld_long = q.mk_lat_data_long([["b", 5]])
assert ld_long.ndim() == 1
assert not ld_long.is_complex()
q.json_results_append(
    f"mk_lat_data_long: ndim={ld_long.ndim()} is_complex={ld_long.is_complex()}"
)

q.json_results_append("set_info with dim_indices")

ld2 = q.LatData()
ld2.set_info([["t", 8], ["mom", 3, ["(0,0,0)", "(1,0,0)", "(1,1,0)"]]])
assert ld2.ndim() == 2
assert ld2.dim_name(1) == "mom"
assert ld2.dim_indices(1) == ["(0,0,0)", "(1,0,0)", "(1,1,0)"]
q.json_results_append(
    f"set_info: ndim={ld2.ndim()} dim_name(1)={ld2.dim_name(1)} indices={ld2.dim_indices(1)}"
)

q.json_results_append("dimension metadata")

ld3 = q.mk_lat_data([["t", 8], ["mom", 3, ["p0", "p1", "p2"]]], is_complex=True)
assert ld3.ndim() == 2
assert ld3.dim_name(0) == "t"
assert ld3.dim_name(1) == "mom"
assert ld3.dim_size(0) == 8
assert ld3.dim_size(1) == 3
assert ld3.dim_names() == ["t", "mom"]
assert ld3.dim_sizes() == [8, 3]
assert ld3.dim_indices(1) == ["p0", "p1", "p2"]
assert ld3.dim_idx(1, "p1") == 1
assert ld3.is_complex()
q.json_results_append(
    f"metadata: ndim={ld3.ndim()} dim_names={ld3.dim_names()} dim_sizes={ld3.dim_sizes()}"
)
q.json_results_append(
    f"metadata: dim_name(0)={ld3.dim_name(0)} dim_size(0)={ld3.dim_size(0)}"
)
q.json_results_append(
    f"metadata: dim_indices(1)={ld3.dim_indices(1)} dim_idx(1,'p1')={ld3.dim_idx(1, 'p1')}"
)
q.json_results_append(f"metadata: is_complex={ld3.is_complex()}")

info_all = ld3.info()
assert len(info_all) == 2
assert info_all[0][0] == "t"
assert info_all[0][1] == 8
q.json_results_append(f"info(): len={len(info_all)} dim0_name={info_all[0][0]}")

info_dim1 = ld3.info(1)
assert info_dim1[0] == "mom"
assert info_dim1[1] == 3
q.json_results_append(f"info(1): name={info_dim1[0]} size={info_dim1[1]}")

ld3.set_dim_name(0, "tau", ["t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7"])
assert ld3.dim_name(0) == "tau"
assert ld3.dim_indices(0) == ["t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7"]
q.json_results_append(f"set_dim_name: dim_name(0)={ld3.dim_name(0)}")

q.json_results_append("data access")

ld4 = q.mk_lat_data([["t", 4], ["x", 3]], is_complex=True)
ld4.set_zero()
ld4[0, 1] = 3.14
val = ld4[0, 1]
assert val == 3.14
q.json_results_append(f"get/set item: ld4[0,1]={val}")

ld4[1, 2] = 1.0 + 2.0j
val_c = ld4[1, 2]
assert val_c == 1.0 + 2.0j
q.json_results_append(f"complex set/get: ld4[1,2]={val_c}")

q.json_results_append("to_numpy / from_numpy")

arr = ld4.to_numpy()
assert arr.shape == (4, 3)
assert arr.dtype == np.complex128
q.json_results_append(f"to_numpy: shape={arr.shape} dtype={arr.dtype}")

ld5 = q.LatData()
ld5.from_numpy(arr)
assert ld5.ndim() == 2
assert np.allclose(np.asarray(ld5), arr)
q.json_results_append(
    f"from_numpy: ndim={ld5.ndim()} match={np.allclose(np.asarray(ld5), arr)}"
)

arr_f = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
ld5b = q.LatData()
ld5b.from_numpy(arr_f, dim_names=["x", "y"], is_complex=False)
assert ld5b.ndim() == 2
assert ld5b.dim_names() == ["x", "y"]
assert not ld5b.is_complex()
q.json_results_append(
    f"from_numpy(dim_names): ndim={ld5b.ndim()} names={ld5b.dim_names()}"
)

q.json_results_append("to_list / from_list")

lst = ld4.to_list()
assert len(lst) == 4 * 3
q.json_results_append(f"to_list: len={len(lst)}")

ld6 = q.LatData()
ld6.from_list(lst)
assert ld6.ndim() == 1
assert ld6.dim_sizes() == [12]
assert np.allclose(np.asarray(ld6).ravel(), np.asarray(ld4).ravel())
q.json_results_append(
    f"from_list: ndim={ld6.ndim()} data_match={np.allclose(np.asarray(ld6).ravel(), np.asarray(ld4).ravel())}"
)

q.json_results_append("save_str / load_str round-trip")

ld7 = q.mk_lat_data([["t", 4], ["p", 2]])
ld7[0, 0] = 1.5
ld7[1, 1] = -2.5
content = ld7.save_str()
assert isinstance(content, bytes)
ld8 = q.LatData()
ld8.load_str(content)
assert np.allclose(np.asarray(ld8), np.asarray(ld7))
q.json_results_append(
    f"save_str/load_str: type={type(content).__name__} match={np.allclose(np.asarray(ld8), np.asarray(ld7))}"
)

q.json_results_append("arithmetic (LatData)")

ld_a = q.mk_lat_data([["t", 4], ["x", 2]])
ld_a[0, 0] = 1.0
ld_a[1, 1] = 2.0

ld_b = q.mk_lat_data([["t", 4], ["x", 2]])
ld_b[0, 0] = 0.5
ld_b[1, 1] = 1.5

ld_sum = ld_a + ld_b
assert ld_sum[0, 0] == 1.5
assert ld_sum[1, 1] == 3.5
q.json_results_append(f"add: (ld_a+ld_b)[0,0]={ld_sum[0, 0]} [1,1]={ld_sum[1, 1]}")

ld_diff = ld_a - ld_b
assert ld_diff[0, 0] == 0.5
assert ld_diff[1, 1] == 0.5
q.json_results_append(f"sub: (ld_a-ld_b)[0,0]={ld_diff[0, 0]} [1,1]={ld_diff[1, 1]}")

ld_scaled = ld_a * 3.0
assert ld_scaled[0, 0] == 3.0
assert ld_scaled[1, 1] == 6.0
q.json_results_append(
    f"mul scalar: (ld_a*3)[0,0]={ld_scaled[0, 0]} [1,1]={ld_scaled[1, 1]}"
)

ld_neg = -ld_a
assert ld_neg[0, 0] == -1.0
assert ld_neg[1, 1] == -2.0
q.json_results_append(f"neg: (-ld_a)[0,0]={ld_neg[0, 0]} [1,1]={ld_neg[1, 1]}")

ld_c = ld_a.copy()
ld_c += ld_b
assert ld_c[0, 0] == 1.5
q.json_results_append(f"iadd: ld_c[0,0]={ld_c[0, 0]}")

ld_c -= ld_b
assert ld_c[0, 0] == 1.0
q.json_results_append(f"isub: ld_c[0,0]={ld_c[0, 0]}")

ld_c *= 2.0
assert ld_c[0, 0] == 2.0
q.json_results_append(f"imul: ld_c[0,0]={ld_c[0, 0]}")

qn = ld_a.qnorm()
assert qn == 1.0 + 4.0
q.json_results_append(f"qnorm: {qn}")

ld_z = ld_a.copy()
ld_z.set_zero()
assert ld_z.qnorm() == 0.0
q.json_results_append(f"set_zero: qnorm={ld_z.qnorm()}")

assert ld_a.is_match(ld_b)
q.json_results_append(f"is_match: {ld_a.is_match(ld_b)}")

q.json_results_append("buffer protocol")

ld_buf = q.mk_lat_data([["t", 3], ["x", 2]])
ld_buf[0, 0] = 5.0
arr_view = np.asarray(ld_buf)
assert arr_view[0, 0] == 5.0
arr_view[1, 1] = 7.0
assert ld_buf[1, 1] == 7.0
q.json_results_append(
    f"buffer protocol: view[0,0]={arr_view[0, 0]} ld[1,1]={ld_buf[1, 1]}"
)

q.json_results_append("serialization (pickle)")

ld_pickle = q.mk_lat_data([["t", 4], ["mom", 2, ["p0", "p1"]]])
ld_pickle[0, 0] = 1.0 + 0.5j
ld_pickle[1, 1] = -3.0
data = pickle.dumps(ld_pickle)
ld_restored = pickle.loads(data)
assert ld_restored.ndim() == ld_pickle.ndim()
assert ld_restored.dim_names() == ld_pickle.dim_names()
assert ld_restored.dim_indices(1) == ["p0", "p1"]
assert np.allclose(np.asarray(ld_restored), np.asarray(ld_pickle))
q.json_results_append(
    f"pickle: ndim={ld_restored.ndim()} names={ld_restored.dim_names()} match={np.allclose(np.asarray(ld_restored), np.asarray(ld_pickle))}"
)

q.json_results_append("copy")

ld_orig = q.mk_lat_data([["t", 3]])
ld_orig[0] = 42.0

ld_copy = ld_orig.copy()
assert np.allclose(np.asarray(ld_copy), np.asarray(ld_orig))
ld_copy[0] = 99.0
assert ld_orig[0] == 42.0
q.json_results_append(f"copy: orig[0]={ld_orig[0]} copy_modified[0]={ld_copy[0]}")

ld_shallow = copy.copy(ld_orig)
assert np.allclose(np.asarray(ld_shallow), np.asarray(ld_orig))
q.json_results_append(
    f"__copy__: match={np.allclose(np.asarray(ld_shallow), np.asarray(ld_orig))}"
)

ld_deep = copy.deepcopy(ld_orig)
assert np.allclose(np.asarray(ld_deep), np.asarray(ld_orig))
q.json_results_append(
    f"__deepcopy__: match={np.allclose(np.asarray(ld_deep), np.asarray(ld_orig))}"
)

q.json_results_append("MPI operations")

ld_mpi = q.mk_lat_data([["t", 4]])
ld_mpi.set_zero()
ld_mpi[0] = 1.0
num_node = q.get_num_node()
ld_mpi_sum = ld_mpi.glb_sum()
assert ld_mpi_sum[0] == float(num_node)
q.json_results_append(f"glb_sum: [0]={ld_mpi_sum[0]} num_node={num_node}")

ld_mpi2 = ld_mpi.copy()
ld_mpi2.glb_sum_in_place()
assert ld_mpi2[0] == float(num_node)
q.json_results_append(f"glb_sum_in_place: [0]={ld_mpi2[0]}")

ld_bc = q.mk_lat_data([["t", 3]])
ld_bc[0] = 5.0
ld_bc.bcast()
assert ld_bc[0] == 5.0
q.json_results_append(f"bcast: [0]={ld_bc[0]}")

q.json_results_append("complex arithmetic")

ld_ca = q.mk_lat_data([["t", 2]])
ld_ca[0] = 1.0 + 2.0j
ld_ca[1] = 3.0 + 4.0j

ld_cb = ld_ca * 2.0
assert ld_cb[0] == 2.0 + 4.0j
q.json_results_append(f"complex mul scalar: {ld_cb[0]}")

ld_cc = ld_ca * (1.0 + 1.0j)
expected = (1.0 + 2.0j) * (1.0 + 1.0j)
assert abs(ld_cc[0] - expected) < 1e-14
q.json_results_append(f"complex mul complex: {ld_cc[0]}")

ld_cn = -ld_ca
assert ld_cn[0] == -(1.0 + 2.0j)
q.json_results_append(f"complex neg: {ld_cn[0]}")

q.json_results_append("LatDataRealF round-trip")

ld_rf = q.mk_lat_data_real_f([["t", 4], ["x", 2]])
assert ld_rf.ndim() == 2
assert ld_rf.is_complex()
ld_rf[0, 0] = 1.5
arr_rf = np.asarray(ld_rf)
assert arr_rf.dtype == np.complex64
content_rf = ld_rf.save_str()
ld_rf2 = q.LatDataRealF()
ld_rf2.load_str(content_rf)
assert np.allclose(np.asarray(ld_rf2), np.asarray(ld_rf))
q.json_results_append(
    f"LatDataRealF: dtype={arr_rf.dtype} round_trip={np.allclose(np.asarray(ld_rf2), np.asarray(ld_rf))}"
)

q.check_log_json(__file__, check_eps=1e-14)
q.timer_display()
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
