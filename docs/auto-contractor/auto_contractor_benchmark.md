# `auto_contractor.benchmark` — Matrix Operation Benchmarks

Source: `qlat/auto_contractor/benchmark.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Benchmark Functions](#benchmark-functions)
3. [Main Block](#main-block)
4. [Examples](#examples)

---

## Overview

`benchmark` provides timer-instrumented wrappers for measuring the throughput
of the low-level matrix operations that dominate auto-contractor evaluation
cost.  Each benchmark function calls its target repeatedly in batches of 10
inside a `q.TimerFork` context so that per-call timing and FLOP rates are
recorded by the Qlattice profiling system.

## Benchmark Functions

### `benchmark_function_1(f, arg, ...)`

Benchmark a single-argument function `f(arg)`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `f` | callable | — | Function to benchmark |
| `arg` | any | — | Argument passed to `f` |
| `benchmark_size` | int | 1000 | Number of outer loop iterations |
| `benchmark_num` | int | 10 | Number of benchmark repetitions |
| `total_flops` | int | 0 | FLOP count reported to the timer |

### `benchmark_function_2(f, arg1, arg2, ...)`

Benchmark a two-argument function `f(arg1, arg2)`.  Parameters are identical
to `benchmark_function_1` except for the extra argument.

### `benchmark_function_3(f, arg1, arg2, arg3, ...)`

Benchmark a three-argument function `f(arg1, arg2, arg3)`.  Same parameter
conventions.

## Main Block

When run as a script, the module benchmarks the core spin-colour matrix
operations used by the auto-contractor runtime:

| Operation | Function | FLOP Count |
|---|---|---|
| `sc * sc` | `mat_mul_sc_sc` | 13536 |
| `sc * s` | `mat_mul_sc_s` | 4320 |
| `s * sc` | `mat_mul_s_sc` | 4320 |
| `s * s` | `mat_mul_s_s` | 4320 |
| `tr(sc, sc)` | `mat_sc_sc_trace` | 480 |
| `tr(sc)` | `mat_sc_trace` | 22 |

Random test matrices are created via `make_rand_spin_color_matrix` and
`make_rand_spin_matrix` with a fixed `RngState("seed")` for reproducibility.

## Examples

```python
import qlat as q
import auto_contractor as ac

q.begin_with_mpi([[1, 1, 1, 4]])

from auto_contractor.benchmark import benchmark_function_2
from auto_contractor.eval import make_rand_spin_color_matrix, get_mat, mat_mul_sc_sc

rng = q.RngState("bench-test")
m1 = make_rand_spin_color_matrix(rng)
m2 = make_rand_spin_color_matrix(rng)

benchmark_function_2(
    mat_mul_sc_sc, get_mat(m1), get_mat(m2),
    benchmark_size=100, benchmark_num=5, total_flops=13536,
)

q.end_with_mpi()
```
