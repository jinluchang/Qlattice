#!/usr/bin/env python3

json_results = []
check_eps = 1e-14

def json_results_append(*args):
    json_results.append(args)

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

cache = q.LRUCache(4)
cache.clear()
v = 12 in cache
json_results_append(f"{v}")
cache[12] = "hello"
v = 12 in cache
json_results_append(f"{v}")
v = cache[12]
json_results_append(f"{v}")
for i in range(8):
    cache[i] = f"{i}"
v = len(cache.cache)
json_results_append(f"{v}")
v = cache.get(12)
json_results_append(f"{v}")
v = cache.get(12, "not found")
json_results_append(f"{v}")
v = cache.cache
json_results_append(f"{v}")

q.check_log_json(__file__, json_results, check_eps=check_eps)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
