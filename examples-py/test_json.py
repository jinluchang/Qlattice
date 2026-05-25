#!/usr/bin/env python3

import numpy as np
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

q.json_results_append("test json round-trip")

# complex
original = complex(1.5, -2.5)
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"complex: {original} -> {restored} ; match={original == restored}"
)

# complex128
original = np.complex128(3.14 + 2.71j)
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"complex128: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# complex64
original = np.complex64(1.0 - 0.5j)
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"complex64: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# float32
original = np.float32(3.14159)
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"float32: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# int64
original = np.int64(9223372036854775807)
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"int64: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# int32
original = np.int32(2147483647)
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"int32: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# ndarray
original = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float64)
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"ndarray: shape={original.shape} -> {restored.shape} ; type={type(restored).__name__} ; match={np.allclose(original, restored)}"
)

# Coordinate
original = q.Coordinate([1, 2, 3, 4])
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"Coordinate: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# CoordinateD
original = q.CoordinateD([1.1, 2.2, 3.3, 4.4])
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"CoordinateD: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# range
original = range(2, 10, 3)
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(f"range: {original} -> {restored} ; match={original == restored}")

# float (plain Python float - should round-trip via standard JSON)
original = 3.14159265358979
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"float: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# int (plain Python int)
original = 42
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(
    f"int: {original} -> {restored} ; type={type(restored).__name__} ; match={original == restored}"
)

# list with mixed types
original = [np.int64(1), np.float32(2.5), complex(3, 4), q.Coordinate([5, 6, 7, 8])]
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(f"list: len={len(restored)}")
q.json_results_append(
    f"list[0]: type={type(restored[0]).__name__} ; match={restored[0] == original[0]}"
)
q.json_results_append(
    f"list[1]: type={type(restored[1]).__name__} ; match={restored[1] == original[1]}"
)
q.json_results_append(
    f"list[2]: type={type(restored[2]).__name__} ; match={restored[2] == original[2]}"
)
q.json_results_append(
    f"list[3]: type={type(restored[3]).__name__} ; match={restored[3] == original[3]}"
)

# dict with mixed types
original = {
    "int64": np.int64(-1),
    "complex128": np.complex128(1 + 2j),
    "coord": q.Coordinate([0, 0, 0, 0]),
    "coord_d": q.CoordinateD([1.5, 2.5, 3.5, 4.5]),
    "arr": np.array([1.0, 2.0, 3.0]),
}
s = q.json_dumps(original)
restored = q.json_loads(s)
q.json_results_append(f"dict: keys={sorted(restored.keys())}")
q.json_results_append(
    f"dict['int64']: type={type(restored['int64']).__name__} ; match={restored['int64'] == original['int64']}"
)
q.json_results_append(
    f"dict['complex128']: type={type(restored['complex128']).__name__} ; match={restored['complex128'] == original['complex128']}"
)
q.json_results_append(
    f"dict['coord']: type={type(restored['coord']).__name__} ; match={restored['coord'] == original['coord']}"
)
q.json_results_append(
    f"dict['coord_d']: type={type(restored['coord_d']).__name__} ; match={restored['coord_d'] == original['coord_d']}"
)
q.json_results_append(
    f"dict['arr']: type={type(restored['arr']).__name__} ; match={np.allclose(restored['arr'], original['arr'])}"
)

q.check_log_json(__file__, check_eps=1e-14)
q.timer_display()
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
