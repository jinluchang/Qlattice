# `qlat_utils.parallel` — Multiprocessing Helpers for Parallel Workloads

Source: `qlat-utils/qlat_utils/parallel.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Core Functions](#core-functions)
3. [Configuration](#configuration)
4. [Utility Functions](#utility-functions)
5. [Examples](#examples)

---

## Overview

`parallel` provides multiprocessing utilities for embarrassingly parallel
workloads in lattice-QCD computations.  It wraps Python's
`multiprocessing.Pool` with:

- Automatic garbage-collection management (freeze/unfreeze) to reduce
  memory overhead in worker processes.
- Process-local initialisation that resets timers and suppresses verbose
  output in workers.
- A lightweight iteration tracer for monitoring long-running map operations.
- A `sum_list` helper for element-wise reduction of list results.

The number of worker processes is controlled by the `q_num_mp_processes`
environment variable (falling back to `q_num_threads` or `OMP_NUM_THREADS`,
default `2`).  Setting it to `0` disables multiprocessing and runs
sequentially.

---

## Core Functions

### `parallel_map`

```python
q.parallel_map(func, iterable, *, n_proc=None, chunksize=1,
               process_initialization=process_initialization, verbose=None)
```

Apply `func` to every element of `iterable` across a pool of `n_proc`
worker processes.

| Parameter | Description |
|---|---|
| `func` | Callable applied to each element (must be picklable) |
| `iterable` | Input iterable |
| `n_proc` | Number of workers (`None` = use env config, `0` = sequential) |
| `chunksize` | Chunk size for `Pool.imap` |
| `process_initialization` | Callable run once per worker on startup |
| `verbose` | Verbosity level (`None` = use env config) |

**Returns:** `list` of results in input order.

When `n_proc=0`, runs as a simple sequential map with optional progress
logging.

### `parallel_map_sum`

```python
q.parallel_map_sum(func, iterable, *, n_proc=None, sum_function=None,
                   sum_start=None, chunksize=1,
                   process_initialization=process_initialization, verbose=None)
```

Like `parallel_map`, but applies `sum_function` to the results iterator
before returning.  This avoids materialising the full result list when only
an aggregate is needed.

| Parameter | Description |
|---|---|
| `sum_function` | Reduction function (default: built-in `sum`) |
| `sum_start` | Optional initial value for the reduction |

**Returns:** The reduced result.

---

## Configuration

### `get_q_num_mp_processes` / `set_q_num_mp_processes`

```python
q.get_q_num_mp_processes()   # reads from env on first call
q.set_q_num_mp_processes(4)  # override at runtime
```

Get or set the default number of worker processes.  On first call, reads
from environment variables in order: `q_num_mp_processes`,
`q_num_threads`, `OMP_NUM_THREADS` (default `2`).

### `get_q_verbose_parallel_map` / `set_q_verbose_parallel_map`

```python
q.get_q_verbose_parallel_map()   # reads from env on first call
q.set_q_verbose_parallel_map(1)  # override at runtime
```

Get or set the default verbosity for `parallel_map` and
`parallel_map_sum`.  On first call, reads from `q_verbose_parallel_map`
(default `2`).

---

## Utility Functions

### `sum_list`

```python
q.sum_list(res, start=None)
```

Element-wise summation of a list of lists.

```python
res = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
q.sum_list(res)  # [12, 15, 18]
```

With an initial value:

```python
q.sum_list(res, start=[10, 10, 10])  # [22, 25, 28]
```

### `trace_iter`

```python
q.trace_iter(iterable, *, tag=None, step_size=None, max_idx=None, verbose_level=0)
```

Generator that yields elements from `iterable` while periodically logging
progress.  Useful for monitoring long-running parallel map results.

| Parameter | Description |
|---|---|
| `tag` | Prefix for log messages (default: function name) |
| `step_size` | Log every `step_size` elements (default: `1`) |
| `max_idx` | Total count for `N/M` progress display |
| `verbose_level` | Logging verbosity level |

---

## Examples

### Basic parallel map

```python
import qlat_utils as q

def square(x):
    return x * x

results = q.parallel_map(square, range(20), n_proc=4)
print(results)  # [0, 1, 4, 9, ..., 361]
```

### Parallel map with sum reduction

```python
import qlat_utils as q

def partial_sum(start_end):
    start, end = start_end
    return sum(range(start, end))

total = q.parallel_map_sum(
    partial_sum,
    [(0, 1000), (1000, 2000), (2000, 3000)],
    n_proc=3,
)
print(total)  # sum(range(3000))
```

### Monitoring with trace_iter

```python
import qlat_utils as q

def heavy_computation(x):
    import time
    time.sleep(0.01)
    return x ** 2

results = q.parallel_map(heavy_computation, range(100), n_proc=4)
for r in q.trace_iter(results, step_size=10, max_idx=100):
    pass  # progress is logged automatically
```

### Sequential fallback (n_proc=0)

```python
import qlat_utils as q

results = q.parallel_map(str, range(5), n_proc=0)
# Runs as simple map(str, range(5)) — no pool created
```
