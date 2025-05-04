#!/usr/bin/env python3

import sys
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

q.json_results_append(f"test LRUCache")

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

q.json_results_append(f"test cache_call")

block_size = 10

block_size_dict = { "48I": 10, }

@q.cache_call(
    maxsize=4,
)
def func(x):
    q.json_results_append(f"run: func({x}) ; block_size={block_size} ; block_size_dict={block_size_dict}")
    return x, block_size, block_size_dict

for i in range(10):
    q.json_results_append(f"i={i} ; func(i % 4)")

@q.cache_call(
    maxsize=4,
    get_state=lambda: (block_size, block_size_dict,),
)
def func(x):
    q.json_results_append(f"run: func({x}) ; block_size={block_size} ; block_size_dict={block_size_dict}")
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
    get_state=lambda: (block_size, block_size_dict,),
    path="cache/func",
)
def func(x):
    q.json_results_append(f"run: func({x}) ; block_size={block_size} ; block_size_dict={block_size_dict}")
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
    q.json_results_append(f"run: func({x}) ; block_size={block_size} ; block_size_dict={block_size_dict}")
    return x, block_size, block_size_dict

for i in range(10):
    q.json_results_append(f"i={i} ; {func(i % 4)}")

with q.TimerFork():
    q.json_results_append(f"run: in TimerFork")

q.check_log_json(__file__, check_eps=1e-14)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
