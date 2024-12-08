#!/usr/bin/env python3
import gpt as g
import numpy as np
conf=g.default.get_int("--conf", None)
if conf is not None:
    rng = g.random(f"96I-positions-{conf}")
    n_all_points = 2654208 # 1/64 of all points
    L = [96,96,96,192]
    all_positions = np.array([
        [rng.uniform_int(min=0, max=L[i] - 1) for i in range(4)] for j in range(n_all_points)
    ], dtype=np.int32)
    print(all_positions)