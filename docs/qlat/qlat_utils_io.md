# `qlat.utils_io` — Lock Management, Shutdown, and Time-Limit Utilities

Source: `qlat/qlat/utils_io.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Lock Management](#lock-management)
3. [Graceful Shutdown](#graceful-shutdown)
4. [Time and Stop Checks](#time-and-stop-checks)
5. [Examples](#examples)

---

## Overview

The `qlat.utils_io` module provides MPI-aware process coordination primitives
for long-running lattice QCD simulations. It wraps C++ utilities for:

- **Lock management** — acquire and release filesystem locks so that multiple
  MPI processes or jobs can coordinate access to shared resources (e.g.,
  gauge configurations on disk).
- **Graceful shutdown** (`qquit`) — clean all Python and C++ caches, then
  terminate the program in an orderly fashion.
- **Time-limit checking** — poll wall-clock budgets or SLURM job end-times to
  allow orderly exit before the scheduler kills the job.
- **Stop-file checking** — poll for the existence of a sentinel file that
  signals the simulation should halt.

All functions in this module are available under the `qlat` (`q`) namespace
after `import qlat as q`.

---

## Lock Management

### `obtained_lock_history_list`

A module-level Python `list` that records every path for which
`obtain_lock` returned `True`. This provides an audit trail of locks
acquired during the lifetime of the process and is used internally by
`release_lock` and `qquit`.

---

### `obtain_lock(path: str) -> bool`

Try to acquire a filesystem lock at `path`. Decorated with `@q.timer`.

| Parameter | Type | Description |
|---|---|---|
| `path` | `str` | Filesystem path for the lock |
| **Returns** | `bool` | `True` if the lock was successfully acquired |

If the lock is acquired, `path` is appended to `obtained_lock_history_list`.
The underlying C++ implementation (`cc.obtain_lock`) uses filesystem-level
locking (typically `mkdir` atomics) that is safe across MPI ranks.

```python
if q.obtain_lock("/scratch/run_001/lock"):
    # this process owns the lock
    do_work()
    q.release_lock()
```

---

### `release_lock()`

Release the currently held lock. Decorated with `@q.timer`.

Calls the C++ `cc.release_lock()` which removes the lock file/directory
created by `obtain_lock`.

---

## Graceful Shutdown

### `qquit(msg: str)`

Clean all Python-level caches (via `q.clean_cache()`), then call the C++
`cc.qquit(msg)` which clears all C++-level caches and terminates the program.

| Parameter | Type | Description |
|---|---|---|
| `msg` | `str` | Message printed before termination |

This is the recommended way to exit an qlat program when an unrecoverable
error is detected or a time/stop condition is met.

```python
q.qquit("finished all trajectories")
```

---

## Time and Stop Checks

<!-- TODO: The return type is not bool. The C++ cc.check_time_limit() returns
     None (falsy) when time is available, or calls qquit() to terminate the
     program when the limit is reached. It never returns True. The doc
     signature and return-type row should reflect this (e.g. "-> None" or
     clarify that it terminates on limit). -->
### `check_time_limit(budget: float = None) -> bool`

Check whether the simulation is approaching its time limit. Decorated with
`@q.timer`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `budget` | `float` or `None` | `None` | Time budget in seconds. If `None`, read from `q.get_time_budget()` |
| **Returns** | `bool` | | `True` if the time limit is approaching |

When `budget` is `None`, the value is read from environment variables (in
order of priority):

| Variable | Description |
|---|---|
| `q_budget` | Budget in seconds, e.g. `export q_budget="$((1 * 60 * 60))"` for 1 hour |
| `q_end_time` | Absolute Unix timestamp, e.g. `export q_end_time="$(($(date +%s) + 12 * 60 * 60))"` or `export q_end_time="$SLURM_JOB_END_TIME"` |

```python
# Check with explicit 30-minute budget
if q.check_time_limit(30 * 60):
    q.qquit("time limit reached")

# Check using environment variable
if q.check_time_limit():
    q.qquit("time limit reached")
```

---

<!-- TODO: The return type is not bool. The C++ cc.check_stop() returns
     None (falsy) when the file does not exist, or calls qquit() to terminate
     the program when the file is found. It never returns True. The doc
     signature and return-type row should reflect this. -->
### `check_stop(fn: str = "stop.txt") -> bool`

Check whether a sentinel file `fn` exists in the current working directory.
Decorated with `@q.timer`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fn` | `str` | `"stop.txt"` | Filename to check |
| **Returns** | `bool` | | `True` if the stop file exists |

This provides a simple mechanism for an operator or batch script to signal a
running simulation to stop gracefully by creating a file.

```python
# In the main simulation loop
for traj in range(start_traj, max_traj):
    do_trajectory(traj)
    if q.check_stop():
        q.qquit("stop file detected")
    if q.check_time_limit():
        q.qquit("time limit reached")
```

---

## Examples

### Lock-Based Coordination

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

path = "/scratch/lattice_run/lock"
if q.obtain_lock(path):
    q.displayln_info("Lock acquired, performing critical section.")
    # ... critical section: read/write shared gauge config ...
    q.release_lock()
    q.displayln_info("Lock released.")
else:
    q.displayln_info("Could not acquire lock, skipping.")

q.end_with_mpi()
```

### Main Loop with Time and Stop Checks

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

max_traj = 1000
for traj in range(max_traj):
    q.displayln_info(f"Starting trajectory {traj}")
    # ... run trajectory ...

    if q.check_stop():
        q.qquit(f"stop file detected at trajectory {traj}")
    if q.check_time_limit():
        q.qquit(f"time limit reached at trajectory {traj}")

q.displayln_info("CHECK: finished successfully.")
q.end_with_mpi()
```

### SLURM-Aware Time Management

```bash
# In your SLURM submission script:
export q_end_time="$SLURM_JOB_END_TIME"
export q_budget="$((30 * 60))"  # 30-minute warning margin

srun python3 run_simulation.py
```

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

# check_time_limit() automatically reads q_end_time and q_budget
# from the environment
if q.check_time_limit():
    q.qquit("SLURM job ending soon, saving and exiting.")

q.end_with_mpi()
```
