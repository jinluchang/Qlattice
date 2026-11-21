#!/usr/bin/env python3

import numpy as np

import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4],
]

q.begin_with_mpi(size_node_list)

q.json_results_append("test LRUCache")

cache = q.LRUCache(4)
cache.clear()
v = 12 in cache
q.json_results_append(f"{v}")
cache[12] = "hello"
v = 12 in cache
q.json_results_append(f"{v}")
v = cache[12]
q.json_results_append(f"{v}")
for i in range(8):
    cache[i] = f"{i}"
v = len(cache.cache)
q.json_results_append(f"{v}")
v = cache.get(12)
q.json_results_append(f"{v}")
v = cache.get(12, "not found")
q.json_results_append(f"{v}")
v = list(cache.cache.items())
q.json_results_append(f"{v}")

q.json_results_append("test cache_call")

block_size = 10

block_size_dict = {
    "48I": 10,
}

@q.cache_call(
    maxsize=4,
)
def func(x):
    q.json_results_append(
        f"run: func({x}) ; block_size={block_size} ; block_size_dict={block_size_dict}"
    )
    return x, block_size, block_size_dict

for i in range(10):
    q.json_results_append(f"i={i} ; func(i % 4)")

@q.cache_call(
    maxsize=4,
    get_state=lambda: (
        block_size,
        block_size_dict,
    ),
)
def func(x):
    q.json_results_append(
        f"run: func({x}) ; block_size={block_size} ; block_size_dict={block_size_dict}"
    )
    return x, block_size, block_size_dict

for i in range(10):
    q.json_results_append(f"i={i} ; func(i % 4)")

block_size = 2

for i in range(10):
    q.json_results_append(f"i={i} ; func(i % 4)")

block_size_dict["64I"] = 20

for i in range(10):
    q.json_results_append(f"i={i} ; func(i % 4)")

for i in range(20):
    q.json_results_append(f"i={i} ; func(i % 10)")

@q.cache_call(
    maxsize=2,
    get_state=lambda: (
        block_size,
        block_size_dict,
    ),
    path="cache/func",
)
def func(x):
    q.json_results_append(
        f"run: func({x}) ; block_size={block_size} ; block_size_dict={block_size_dict}"
    )
    return x, block_size, block_size_dict

for i in range(10):
    q.json_results_append(f"i={i} ; func(i % 4)")

block_size = 2

for i in range(10):
    q.json_results_append(f"i={i} ; func(i % 4)")

block_size_dict["64I"] = 20

for i in range(20):
    q.json_results_append(f"i={i} ; func(i % 10)")

@q.cache_call(
    maxsize=4,
    get_state=lambda: (block_size,),
    is_hash_args=False,
)
def func(x):
    q.json_results_append(
        f"run: func({x}) ; block_size={block_size} ; block_size_dict={block_size_dict}"
    )
    return x, block_size, block_size_dict

for i in range(10):
    q.json_results_append(f"i={i} ; {func(i % 4)}")

with q.TimerFork():
    q.json_results_append("run: in TimerFork")

# Regression asserts for matrix scalar arithmetic / reflected ops / g5_herm.
# Asserts only: nothing is printed or recorded, so the reference logs are unchanged.
wm = q.WilsonMatrix()
wm[:] = 0.0
wm[0, 0] = 1.0 + 2.0j
wm[1, 0] = 3.0 - 4.0j
wm_ref = wm.copy()
assert 2 * wm == wm * 2
assert np.float64(2.0) * wm == wm * 2
assert 2 * wm == q.mat_mul_a_wm(2.0, wm)
assert wm == wm_ref
sm = q.SpinMatrix()
sm[:] = 0.0
sm[0, 0] = 1.0 + 2.0j
sm_ref = sm.copy()
assert 2 * sm == sm * 2
assert sm == sm_ref
cm = q.ColorMatrix()
cm[:] = 0.0
cm[0, 0] = 1.0 + 2.0j
cm_ref = cm.copy()
assert 2 * cm == cm * 2
assert cm == cm_ref
wm_g5_ref = q.wilson_matrix_g5_herm(wm_ref)
assert wm.g5_herm() is wm
assert wm == wm_g5_ref
for m in (wm, sm, cm):
    try:
        hash(m)
    except TypeError:
        pass
    else:
        raise AssertionError("matrix types must be unhashable")

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-14)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
