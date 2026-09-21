# `qlat_utils.jackknife_utils` — Jackknife Resampling

Source: `qlat-utils/qlat_utils/jackknife_utils.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Jackknife Resampling](#jackknife-resampling)
3. [Super-Jackknife](#super-jackknife)
4. [Randomized Jackknife-Bootstrap Hybrid](#randomized-jackknife-bootstrap-hybrid)
5. [Unified Jackknife API](#unified-jackknife-api)
6. [Context Managers](#context-managers)
7. [Examples](#examples)

---

## Overview

The `qlat_utils.jackknife_utils` module provides jackknife resampling and error
estimation for lattice QCD measurements. It includes:

- **Jackknife resampling** — the plain jackknife of a data list, and the
  `jk_avg` / `jk_err` extraction from the resulting array.
- **Super-jackknife** — resampling over the jackknife blocks of the data, with
  the `sjk_*` extraction functions.
- **Randomized jackknife-bootstrap hybrid (RJK)** — random samples with a
  blocking dependent shift, the default of the unified API, with the `rjk_*`
  extraction functions.
- **The unified `g_*` API** — `g_mk_jk` and friends, configured by
  `default_g_jk_kwargs` (or the `JkKwargs` context manager), which dispatch on
  `jk_type` and cache their intermediate state through `get_jk_state` /
  `set_jk_state`.
- **Collective MPI variants** — `is_sync_node=True` (every node holds the whole
  input and obtains the whole result) and `g_mk_jk_distributed` (the input is
  split between the nodes), built on the `mpi4py` collectives.

The generic helpers used here (`q`, the type tuples, `use_kwargs`, `average`,
`block_data`, `filter_np_results`, `fsqr`, `fsqrt`, `qnorm` and
`NewDictValues`) live in `qlat_utils.data`.

```python
import qlat_utils as q

data = [1.0, 1.1, 0.9, 1.05]
jk_arr = q.g_mk_jk(data, list(range(len(data))))
avg, err = q.g_jk_avg_err(jk_arr)
```

---

## Jackknife Resampling

### `jackknife(data_list, *, eps=1)`

Perform standard jackknife. Returns `jk_arr` of length `N + 1` where:
- `jk_arr[0]` = average
- `jk_arr[i]` = `avg - (eps / N) * (data[i] - avg)` for `i >= 1`

### `jk_avg(jk_arr)`

Return the average (first element) of a jackknife array.

### `jk_err(jk_arr, *, eps=1, block_size=1)`

Return the jackknife error estimate:

$$\frac{1}{\text{eps}} \sqrt{ \frac{N}{N - \text{block\_size}} \sum_{i=1}^{N} (jk[i] - \text{jk\_avg})^2 }$$

The `eps` and `block_size` must match those used in the corresponding
`jackknife` call. Note: `len(jk_arr) = N + 1`.

### `jk_avg_err(jk_arr, *, eps=1, block_size=1)`

Return `(jk_avg, jk_err)`.

---

## Super-Jackknife

### `sjackknife(data_list, jk_idx_list, *, avg=None, ...)`

Perform super-jackknife resampling. Data from different ensembles (identified
by `jk_idx_list`) are combined into a single jackknife array.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data_list` | list/ndarray | — | Original data |
| `jk_idx_list` | list | — | Index for each data point (e.g., `(job_tag, traj)`) |
| `avg` | any | `None` | Pre-computed average (auto-computed if `None`) |
| `is_hash_jk_idx` | bool | `True` | Use hash when jk_idx not in `all_jk_idx` |
| `jk_idx_hash_size` | int | `1024` | Hash table size |
| `rng_state` | `RngState` | `None` | RNG state (default: `RngState("rejk")`) |
| `all_jk_idx` | list | `None` | All possible indices; `all_jk_idx[0]` must be `"avg"` |
| `get_all_jk_idx` | callable | `None` | Function returning `all_jk_idx` |
| `jk_blocking_func` | callable | `None` | `(i, jk_idx) -> blocked_jk_idx` |
| `eps` | float | `1` | Scaling factor |

### `sjackknife_distributed(data_list, jk_idx_list, *, avg=None, ...)`

`sjackknife` for a data set that is split between the MPI nodes: `data_list`
and `jk_idx_list` are the local parts of the data set and the returned
`jk_arr` is the local part of the result, so that

```python
jk_arr = np.concatenate(q.get_comm().allgather(jk_local))
```

is the complete data set (node `r` owns the samples in
`range(*get_distributed_range(len(all_jk_idx), r, num_node))`). It is a
collective operation and every node must call it with the same parameters.
The result agrees with `sjackknife` up to the floating-point roundoff but not
bit-for-bit. `g_mk_jk_distributed` dispatches to it.

### `sjackknife_sync_node(data_list, jk_idx_list, *, avg=None, ...)`

`sjackknife` as a collective MPI operation where every node has the whole
input: the input is split between the nodes, `sjackknife_distributed` is
called on the local parts and the parts are gathered, so that every node
obtains the complete `jk_arr`. The result agrees with `sjackknife` up to the
floating-point roundoff but not bit-for-bit. `g_mk_jk_sync_node` dispatches to
it.

### `sjk_avg(jk_arr)` / `sjk_err(jk_arr, *, eps=1)` / `sjk_avg_err(jk_arr, *, eps=1)`

Average, error, and `(avg, err)` for super-jackknife arrays. The error formula
differs from standard jackknife: no `N/(N-1)` factor.

### `sjk_mk_jk_val(rs_tag, val, err, *, ...)`

Create a synthetic jackknife array from a central value and error using
Gaussian random numbers.

---

## Randomized Jackknife-Bootstrap Hybrid

### `rjackknife(data_list, jk_idx_list, *, avg=None, ...)`

Jackknife-bootstrap hybrid resampling. Returns `jk_arr` of length
`1 + n_rand_sample`. The distribution of `jk_arr` approximates the
distribution of the average.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data_list` | list/ndarray | — | Original data |
| `jk_idx_list` | list | — | Index for each data point |
| `avg` | any | `None` | Pre-computed average |
| `rng_state` | `RngState` | `None` | RNG state |
| `n_rand_sample` | int | `1024` | Number of random samples |
| `jk_blocking_func` | callable | `None` | `(i, jk_idx) -> blocked_jk_idx` |
| `is_normalizing_rand_sample` | bool | `False` | Normalize random vectors |
| `is_apply_rand_sample_jk_idx_blocking_shift` | bool | `True` | Shift blocking per sample |
| `eps` | float | `1` | Scaling factor |
| `is_sync_node` | bool | `False` | Collective MPI operation (see below) |

The formula is:

$$jk\_arr[i] = \text{avg} + \sum_{j=1}^{N} \frac{-\text{eps}}{\sqrt{N(N - b(i,j))}} r_{i,j} (d_j - \text{avg})$$

where $r_{i,j} \sim \mathcal{N}(0, 1)$ and $b(i,j)$ is the block size.

If `is_sync_node` is `True`, the operation is assumed to be a collective
operation in a MPI program where every node has the same input. It is then
performed by `rjackknife_sync_node`: the input is split between the nodes and
`rjackknife_distributed` is called on the local parts; the parts are gathered
with `mpi4py`'s `Allgatherv` on the communicator returned by `q.get_comm()`, so
that every node obtains the complete `jk_arr`. The result agrees with the
`is_sync_node=False` result up to the floating-point roundoff, but not
bit-for-bit, because the average and the sums over the data set are reduced
across the nodes. The input is split using the qlat node numbering
(`q.get_id_node()`), which for some Grid processor layouts differs from the
`MPI_COMM_WORLD` rank — `q.get_comm()` is the communicator that matches it (it
must satisfy `comm.rank == q.get_id_node()`). `qlat` and `mpi4py` are imported
only when `is_sync_node` is `True`. This requires qlat to be initialized on the
whole MPI communicator (`q.begin_with_mpi()`, `q.begin_with_gpt()` or
`q.begin_with_grid()`), and the data must have a numeric dtype supported by
MPI.

### `rjackknife_distributed(data_list, jk_idx_list, *, avg=None, ...)`

`rjackknife` for a data set that is split between the MPI nodes: `data_list`
and `jk_idx_list` are the local parts of the data set and the returned
`jk_arr` is the local part of the result, so that

```python
jk_arr = np.concatenate(q.get_comm().allgather(jk_local))
```

is the complete data set (node `r` owns the samples in
`range(*get_distributed_range(1 + n_rand_sample, r, num_node))`). It is a
collective operation and every node must call it with the same parameters.
The result agrees with `rjackknife` up to the floating-point roundoff but not
bit-for-bit. `g_mk_jk_distributed` dispatches to it.

### `rjackknife_sync_node(data_list, jk_idx_list, *, avg=None, ...)`

`rjackknife` as a collective MPI operation where every node has the whole
input: the input is split between the nodes, `rjackknife_distributed` is
called on the local parts and the parts are gathered, so that every node
obtains the complete `jk_arr`. The result agrees with `rjackknife` up to the
floating-point roundoff but not bit-for-bit. `g_mk_jk_sync_node` dispatches to
it.

### `rjk_avg(jk_arr)` / `rjk_err(jk_arr, eps=1)` / `rjk_avg_err(rjk_list, eps=1)`

Average, error, and `(avg, err)` for randomized jackknife arrays.

### `rjk_mk_jk_val(rs_tag, val, err, *, ...)`

Create a synthetic RJK array from a central value and error.

---

## Unified Jackknife API

The `g_*` functions provide a unified interface that dispatches to either
super-jackknife or RJK based on global settings in `default_g_jk_kwargs`.

### `default_g_jk_kwargs`

Global dictionary controlling jackknife behavior. Key settings:

| Key | Default | Description |
|---|---|---|
| `jk_type` | `"rjk"` | `"rjk"` or `"super"` |
| `eps` | `1` | Scaling factor |
| `n_rand_sample` | `1024` | Number of random samples (`"rjk"` only) |
| `is_normalizing_rand_sample` | `False` | Normalize the random vectors (`"rjk"` only) |
| `is_apply_rand_sample_jk_idx_blocking_shift` | `True` | Shift the blocking for every random sample (`"rjk"` only) |
| `is_use_old_rand_alg` | `False` | Reproduce the old random numbers (`"v1"`, `"rjk"` only) |
| `is_hash_jk_idx` | `True` | Hash the `jk_idx` which are not in `all_jk_idx` (`"super"` only) |
| `jk_idx_hash_size` | `1024` | Number of hash based samples (`"super"` only) |
| `all_jk_idx` | `None` | The samples of `"super"`; `None` means the `jk_idx_hash_size` hash samples |
| `get_all_jk_idx` | `None` | Callable returning `all_jk_idx` |
| `block_size` | `1` | Default blocking size |
| `block_size_dict` | `{"job_tag": 1}` | Per-`job_tag` blocking sizes |
| `jk_blocking_func` | `jk_blocking_func_default` | `(i, jk_idx) -> blocked_jk_idx` |
| `rng_state` | `RngState("rejk")` | RNG state (random numbers of `"rjk"`) |
| `is_sync_node` | `False` | Run as a collective MPI operation |

`is_sync_node` is not part of `get_jk_state` / `set_jk_state`: it changes only
how the result is computed (the result agrees up to the floating-point
roundoff), so it is deliberately excluded from the `q.cache_call` cache key.

### `g_mk_jk(data_list, jk_idx_list, *, avg=None, ...)`

Create a (randomized) super-jackknife data set from un-jackknifed data:
`jk_arr[0]` is the average of the data and `jk_arr[1:]` are the resampled
values used by `g_jk_avg_err`. The data set has `g_jk_size()` samples
(`1 + n_rand_sample` for `"rjk"` and `1 + len(all_jk_idx)` for `"super"`) and
the dtype of the data. It dispatches to `sjackknife` or `rjackknife` based on
`jk_type`, and to `g_mk_jk_sync_node` or `g_mk_jk_distributed` for the
collective MPI modes.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data_list` | list/ndarray | — | The un-jackknifed data (`float`, `complex` or `np.ndarray` each); `None` entries are ignored |
| `jk_idx_list` | list | — | Names of the entries of `data_list` (e.g. `(job_tag, traj)`), same length as `data_list` |
| `avg` | any | `None` | Average of the whole data set; computed from `data_list` when `None` |
| `is_sync_node` | bool | `False` | Collective MPI operation with the whole input on every node (see `g_mk_jk_sync_node`) |

All other keyword parameters are the shared settings listed under
`default_g_jk_kwargs` above, e.g. `jk_type`, `eps`, `n_rand_sample`,
`block_size`/`block_size_dict`, `all_jk_idx` and `rng_state`. When the
`data_list` is already jackknifed, increase `eps` by the factor
`len(data_list)`.

### `g_mk_jk_sync_node(data_list, jk_idx_list, *, avg=None, ...)`

Create a (randomized) super-jackknife data set as a collective MPI operation
in which every node holds the whole input and obtains the whole result; this is
what `g_mk_jk(..., is_sync_node=True)` calls. Every node must call it with the
same parameters and with the complete `data_list`/`jk_idx_list`, and every node
obtains the complete `jk_arr` in the same order as the one from `g_mk_jk` on a
single node. The input is split between the nodes in the order of the nodes
(`get_distributed_range`), the local parts are computed by
`rjackknife_sync_node` or `sjackknife_sync_node` and the parts are gathered, so
the work is parallelized over the nodes.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data_list` | list/ndarray | — | The whole data set, the same on every node |
| `jk_idx_list` | list | — | Indices of the whole data set, the same on every node |
| `avg` | any | `None` | Average of the whole data set, the same on every node; computed with `Allreduce` when `None` |

Returns the complete `jk_arr` on every node. Both `jk_type == "rjk"` and
`jk_type == "super"` are supported, and the remaining keyword parameters are
the same as for `g_mk_jk`. Requires the qlat communicator initialized on the
whole MPI communicator (`q.begin_with_mpi()`, `q.begin_with_gpt()`,
`q.begin_with_grid()` or `q.set_comm(...)`) and an MPI-supported data dtype.
The result agrees with `g_mk_jk` up to the floating-point roundoff, but not
bit-for-bit.

### `g_mk_jk_distributed(data_list, jk_idx_list, *, avg=None, ...)`

Create a (randomized) super-jackknife data set when the data set itself is
split between the MPI nodes. This is a collective MPI operation: every node
must call it with the same parameters and with its own disjoint part of the
data set (`data_list` and `jk_idx_list` are the local parts, and together they
must cover the whole data set exactly once), and every node returns its own
part of the result:

```python
import numpy as np
import qlat as q

comm = q.get_comm()
i_start, i_end = q.get_distributed_range(len(data_list), comm.rank, comm.size)
jk_local = q.g_mk_jk_distributed(
    data_list[i_start:i_end], jk_idx_list[i_start:i_end]
)
jk_arr = np.concatenate(comm.allgather(jk_local))  # the complete data set
```

The split of the input between the nodes is free as long as the local parts are
disjoint and cover the whole data set; the split of the output is fixed: the
samples are split between the nodes in the order of the nodes, i.e. node `r`
out of `num_node` owns the samples in
`range(*get_distributed_range(g_jk_size(), r, num_node))` (a node owns 0 samples
when there are more nodes than samples), so concatenating the local parts in
the order of the nodes reproduces the complete data set in the same order as
the one obtained with `g_mk_jk`. It dispatches to `rjackknife_distributed` for
`jk_type == "rjk"` and to `sjackknife_distributed` for `jk_type == "super"`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data_list` | list/ndarray | — | The local part of the un-jackknifed data |
| `jk_idx_list` | list | — | Indices of the local part, same length as the local `data_list` |
| `avg` | any | `None` | Average of the whole data set, the same on every node; computed with `Allreduce` when `None` |

Returns the local part of `jk_arr` (an array with 0 samples on a node which
owns no sample). The remaining keyword parameters are the same as for
`g_mk_jk`. Requires the qlat communicator initialized on the whole MPI
communicator and an MPI-supported data dtype. Only the small `jk_idx` metadata
is gathered on every node, so no node needs to hold the whole data set, the
whole random matrix or the whole result.

The result agrees with `g_mk_jk` (and with `g_mk_jk(..., is_sync_node=True)`)
up to the floating-point roundoff, but not bit-for-bit: the average and the
sums over the data set are reduced across the nodes, which changes the order of
the floating-point additions. The `is_sync_node` setting of
`default_g_jk_kwargs` is ignored, since the output is always distributed.

### `get_distributed_range(total_size, id_node, num_node)`

Return `(i_start, i_end)`, the contiguous range of the `total_size` indices
owned by node `id_node` out of `num_node` nodes:

```
i_start = (total_size * id_node) // num_node
i_end   = (total_size * (id_node + 1)) // num_node
```

The parts are ordered by `id_node`, cover `range(total_size)` exactly once and
differ in size by at most one.

### `get_collective_comm(tag)` / `is_mpi_dtype(dtype)`

Helpers for the collective MPI operations. `get_collective_comm(tag)` returns
`(comm, id_node, num_node)` from `q.get_comm()` and raises an `Exception`
(with `tag` naming the caller) when the qlat communicator is missing or does
not cover the whole MPI communicator; it asserts `comm.rank ==
q.get_id_node()`. `is_mpi_dtype(dtype)` reports whether a numpy dtype can be
summed by MPI.

### `get_distributed_jk_input(fname, comm, data_list, jk_idx_list)`

Return the common state of the distributed jackknife functions:
`(data_arr, jk_idx_list_local, jk_idx_list_glb, n, elem_shape, dtype)`, i.e.
the local data (shape `(n_local, *elem_shape)`), the local and the gathered
`jk_idx`, the total number of data points and the description of a single data
point; `n` is reduced and `jk_idx_list_glb` is gathered over the nodes.

### `get_distributed_avg(comm, data_arr, n, avg=None)`

Return the average of the whole distributed data set: `avg` itself when it is
not `None`, otherwise the local sums reduced with `Allreduce` and divided by
`n`.

### `get_reduce_scattered_jk_arr(partial_arr, comm, id_node, num_node, avg)`

Return the local part of the result from `partial_arr` (shape
`(total_size, *elem_shape)`, the contribution of the local data to every
sample, with the sample 0 equal to `avg` on its owner): the contributions are
summed with `Reduce_scatter` and `avg` is added to the samples. `partial_arr`
may be any array, in particular a non-contiguous view such as a column of a 2-D
array (a contiguous copy is made when needed, since the buffer of the
collective must be contiguous).

### `get_gathered_jk_arr(jk_local, comm, num_node)`

Return the complete result from the local parts `jk_local`, gathered with
`Allgatherv` in the order of the nodes. `jk_local` may be any array, in
particular a non-contiguous view such as a column of a 2-D array (a contiguous
copy is made when needed, since the send buffer of the collective must be
contiguous).

### `g_mk_jk_val(rs_tag, val, err, *, ...)`

Create a synthetic jackknife array from a value and error. Dispatches to
`sjk_mk_jk_val` or `rjk_mk_jk_val`.

### `g_jk_avg(jk_arr)` / `g_jk_err(jk_arr)` / `g_jk_avg_err(jk_arr)`

Unified average, error, and `(avg, err)` extraction.

### `g_jk_avg_err_arr(jk_arr)`

Return an array with shape `jk_arr[0].shape + (2,)` where the last axis is
`(avg, err)`.

### `g_jk_size(*, jk_type, ...)`

Return the number of samples in the jackknife array (`1 + n_samples`).

### `g_jk_blocking_func(i, jk_idx)`

Apply the configured blocking function.

### `g_jk_sample_size(job_tag, traj_list)`

Return the number of distinct blocks for a given `job_tag` and trajectory list.

### `get_jk_state()` / `set_jk_state(state)`

Save and restore the current `default_g_jk_kwargs` state (for use with
`@cache_call`).


---

## Context Managers

### `class JkKwargs(**kwargs)`

Context manager that temporarily overrides `default_g_jk_kwargs`.

```python
with q.JkKwargs(n_rand_sample=2048, block_size=10):
    jk_arr = q.g_mk_jk(data_list, jk_idx_list)
```

---


## Examples

### Jackknife Resampling

```python
import qlat_utils as q
import numpy as np

data = [1.0, 1.1, 0.9, 1.05, 0.95, 1.02, 0.98, 1.03]

# Standard jackknife
jk_arr = q.jackknife(data)
avg = q.jk_avg(jk_arr)
err = q.jk_err(jk_arr)
print(f"Jackknife: avg = {avg:.4f}, err = {err:.4f}")

# Unified API (uses RJK by default)
jk_arr = q.g_mk_jk(data, list(range(len(data))))
avg, err = q.g_jk_avg_err(jk_arr)
print(f"RJK: avg = {avg:.4f}, err = {err:.4f}")
```

### Distributed MPI

```python
import numpy as np
import qlat_utils as q

# Every node holds the whole data set and obtains the whole result.
jk_arr = q.g_mk_jk(data_list, jk_idx_list, is_sync_node=True)

# Or let every node contribute its own part of the data set and of the result.
comm = q.get_comm()
i_start, i_end = q.get_distributed_range(len(data_list), comm.rank, comm.size)
jk_local = q.g_mk_jk_distributed(
    data_list[i_start:i_end], jk_idx_list[i_start:i_end]
)
jk_arr = np.concatenate(comm.allgather(jk_local))
```
