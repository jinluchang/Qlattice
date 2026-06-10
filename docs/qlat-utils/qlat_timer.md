# `qlat_utils.timer` — Timing, Profiling, and Display Utilities

Source: `qlat-utils/qlat_utils/timer.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Verbose Level Control](#verbose-level-control)
3. [Display Functions](#display-functions)
4. [Timer Class](#timer-class)
5. [TimerNone Class](#timernone-class)
6. [Node and MPI Helpers](#node-and-mpi-helpers)
7. [Time Functions](#time-functions)
8. [Timer Display and Reset Functions](#timer-display-and-reset-functions)
9. [TimerFork Context Manager](#timerfork-context-manager)
10. [Timer Decorators](#timer-decorators)
11. [Examples](#examples)

---

## Overview

`timer` provides the profiling and diagnostic infrastructure for `qlat_utils`:

- **Verbose-level gated printing** — `display`, `displayln`, `display_info`, `displayln_info` only emit output when the current verbose level meets a threshold.
- **`Timer` / `TimerNone`** — C++ backed scope-based timers that accumulate call counts, wall-clock time, and optionally FLOP counts.
- **`@timer` decorator family** — one-line function-level profiling with optional FLOP tracking and `TimerFork` integration.
- **Time management** — access to wall-clock timers, configurable time limits, and remaining-time queries for graceful shutdown.
- **`TimerFork`** — context manager that isolates timer statistics and optionally adjusts verbose level for a code block.

Python access:

```python
import qlat_utils as q
```

---

## Verbose Level Control

### `get_verbose_level() -> int`

Return the current verbosity level. Default depends on the `q_verbose` environment variable; if unset, the default is `-1`.

### `set_verbose_level(level=None)`

Set the verbosity level. If `level` is `None`, reset to the default value.

```python
q.set_verbose_level(0)   # quiet
q.set_verbose_level(10)  # very verbose
q.set_verbose_level()    # reset to default
```

---

## Display Functions

### `display(level, *args)`

Print arguments without a trailing newline. If `level` is an `int`, only print when `level <= get_verbose_level()`. If `level` is not an integer, always print. Prefixes output with `Qlat <elapsed_time>:`.

### `display_info(*args)`

Same as `display`, but only prints on MPI rank 0 (`get_id_node() == 0`).

### `displayln(level, *args)`

Same as `display` but appends a newline.

### `displayln_info(*args)`

Same as `displayln`, but only prints on MPI rank 0. This is the primary logging function used throughout qlat for CHECK lines in tests.

### `set_display_method(method=None)`

Set the underlying C++ display function. `method=None` restores the default; `method="py_stdout"` routes output through Python's `sys.stdout`.

---

## Timer Class

```cython
cdef class Timer:
    def __cinit__(self, const cc.std_string& fname, cc.Bool is_verbose=False)
    def start(self)
    def stop(self)
```

A C++ backed scope timer. Each `Timer` accumulates:

| Metric | Description |
|---|---|
| Call count | Number of `start()`/`stop()` pairs |
| Wall-clock time | Total elapsed time between `start()` and `stop()` |
| FLOP count | Optionally incremented externally via `qtimer.flops += n` |

### Parameters

| Parameter | Type | Description |
|---|---|---|
| `fname` | `str` | Timer label, displayed in timer reports |
| `is_verbose` | `bool` | If `True`, print start/stop messages when verbose level is high enough |

### Methods

| Method | Description |
|---|---|
| `start()` | Begin timing. Optionally prints a start message. |
| `stop()` | End timing and accumulate elapsed time. Optionally prints a stop message. |

---

## TimerNone Class

```python
class TimerNone:
    def start(self):  pass
    def stop(self):   pass
```

A no-op drop-in replacement for `Timer`. Useful when timing is conditionally disabled but the call-site code should remain unchanged.

---

## Node and MPI Helpers

### `get_id_node() -> int`

Return the MPI rank (node id). Always returns `0` if MPI has not been initialized.

### `get_num_node() -> int`

Return the total number of MPI ranks. Always returns `1` if MPI has not been initialized.

### `sync_node(tag=None)`

Perform a global synchronization (via `glb_sum`) to ensure all MPI ranks are synced.

---

## Time Functions

### `get_time() -> float`

Return current wall-clock time in seconds since epoch.

### `get_start_time() -> float`

Return the start time (reset by `timer_reset`).

### `get_actual_start_time() -> float`

Return the actual start time (not reset by `timer_reset`).

### `get_total_time() -> float`

Return elapsed time in seconds since the last `timer_reset`.

### `get_actual_total_time() -> float`

Return elapsed time in seconds since actual program start (not reset by `timer_reset`).

### `get_time_limit() -> float`

Return the configured program time limit in seconds. Used by `get_remaining_time()`.

### `set_time_limit(time_limit=None)`

Set the program time limit in seconds. If `None`, reset to the default value.

### `get_remaining_time() -> float`

Return `time_limit - get_total_time()`, i.e. remaining seconds before the time limit is reached.

### `get_time_budget() -> float` / `set_time_budget(time_budget=None)`

Get or set the time budget in seconds.

---

## Timer Display and Reset Functions

### `timer_display(tag="")`

Display all accumulated timer statistics, tagged with `tag`.

### `timer_autodisplay()`

Display timers that are configured for automatic display.

### `timer_display_stack()` / `timer_display_stack_always()`

Display the current timer call stack (for debugging nested timer usage).

### `timer_reset(max_call_times_for_always_show_info=-1)`

Reset all timers, `get_total_time`, and `get_start_time`. Does **not** reset `get_actual_start_time` or `get_actual_total_time`.

### `timer_fork(max_call_times_for_always_show_info=-1)` / `timer_merge()`

Fork/merge the global timer tree. `timer_fork` pushes a new timer context; `timer_merge` pops it and accumulates statistics into the parent.

---

## TimerFork Context Manager

```python
class TimerFork:
    def __init__(self, *, show_display=10, display_tag="TimerFork",
                 verbose=None, max_call_times_for_always_show_info=-1)
```

A context manager that calls `timer_fork()` on entry and `timer_merge()` on exit. The verbose level is always restored to its original value on exit.

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `show_display` | `int` or `bool` | `10` | If `bool`, always/never display. If `int`, display only when `get_total_time() >= show_display` (seconds). |
| `display_tag` | `str` | `"TimerFork"` | Tag prefix for display output |
| `verbose` | `int` or `None` | `None` | If set, temporarily change verbose level inside the block |
| `max_call_times_for_always_show_info` | `int` | `-1` | Passed to `timer_fork` |

### Usage

```python
with q.TimerFork(verbose=5):
    heavy_computation()
# timers auto-display if elapsed >= 10s
```

---

## Timer Decorators

### `@timer` / `@timer(**kwargs)`

The main function-timing decorator.

### Keyword Arguments

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fname` | `str` or `None` | `None` | Custom timer name. Defaults to `"py:<qualname>"`. |
| `is_verbose` | `bool` | `False` | Print start/stop messages |
| `is_flops` | `bool` | `False` | If `True`, function must return `(flops, ret)`; decorator accumulates flops |
| `is_timer_fork` | `bool` | `False` | If `True`, wrap each call in a `TimerFork` |
| `timer_fork_verbose` | `int` or `None` | `None` | Verbose level for the `TimerFork` |
| `timer_fork_max_call_times_for_always_show_info` | `int` | `-1` | Passed to `TimerFork` |
| `timer_fork_show_display` | `int` | `10` | Seconds threshold for auto-display in `TimerFork` |

```python
@q.timer
def compute(x):
    return x * x

@q.timer(is_verbose=True)
def debug_compute(x):
    return x * x
```

### `@timer_fname(fname)`

Shorthand for `@timer(fname=fname)`.

### `@timer_verbose`

Shorthand for `@timer(is_verbose=True)`.

### `@timer_verbose_fname(fname)`

Shorthand for `@timer(fname=fname, is_verbose=True)`.

### `@timer_flops`

Decorator for functions that return `(flops, ret)`. Accumulates the flops count and returns only `ret`.

### `@timer_flops_fname(fname)`

Shorthand for `@timer(fname=fname, is_flops=True)`.

### `@timer_verbose_flops`

Shorthand for `@timer(is_verbose=True, is_flops=True)`.

### `@timer_verbose_flops_fname(fname)`

Shorthand for `@timer(fname=fname, is_verbose=True, is_flops=True)`.

### `default_timer_kwargs`

A `dict` holding the default values for all `@timer` keyword arguments. Useful for introspection.

---

## Examples

### Basic Function Timing

```python
import qlat_utils as q

@q.timer
def slow_function(n):
    total = 0
    for i in range(n):
        total += i
    return total

slow_function(1000000)
q.timer_display()
```

### Verbose Timing with Custom Name

```python
import qlat_utils as q

@q.timer_fname("matrix-inversion")
def invert(m):
    ...

invert(m)
q.timer_display()
```

### FLOP Tracking

```python
import qlat_utils as q

@q.timer_flops
def matmul(a, b):
    n = a.shape[0]
    result = a @ b
    flops = 2 * n**3
    return flops, result

result = matmul(a, b)  # result is just the product, flops tracked internally
```

### Verbose Level Gated Logging

```python
import qlat_utils as q

q.set_verbose_level(0)  # quiet mode

q.displayln_info(0, "this prints on rank 0")   # prints
q.displayln_info(10, "this is hidden")          # suppressed (10 > 0)
```

### TimerFork for Section Profiling

```python
import qlat_utils as q

with q.TimerFork(verbose=2, show_display=5):
    step1()
    step2()
# If block takes >= 5 seconds, timer stats are displayed automatically
```

### Remaining Time Check for Graceful Shutdown

```python
import qlat_utils as q

q.set_time_limit(3600)  # 1 hour

for i in range(1000):
    if q.get_remaining_time() < 60:
        q.displayln_info(0, "Less than 60s remaining, saving and exiting.")
        save_checkpoint()
        break
    do_iteration(i)
```
