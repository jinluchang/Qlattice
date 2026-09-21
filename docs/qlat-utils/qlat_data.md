# `qlat_utils.data` — Data Analysis Utilities

Source: `qlat-utils/qlat_utils/data.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Physical Constants](#physical-constants)
3. [Type Tuples](#type-tuples)
4. [Interpolation](#interpolation)
5. [Data Wrapper Class](#data-wrapper-class)
6. [Basic Statistics](#basic-statistics)
7. [Jackknife Resampling](#jackknife-resampling)
8. [Super-Jackknife](#super-jackknife)
9. [Randomized Jackknife-Bootstrap Hybrid](#randomized-jackknife-bootstrap-hybrid)
10. [Unified Jackknife API](#unified-jackknife-api)
11. [Value Display](#value-display)
12. [Context Managers](#context-managers)
13. [Examples](#examples)

---

## Overview

The `qlat_utils.data` module provides data analysis tools for lattice QCD
computations. It includes:

- **Interpolation** — linear interpolation for arrays with fractional indices.
- **Jackknife resampling** — standard, super-jackknife, and randomized
  jackknife-bootstrap (RJK) methods for error estimation.
- **The `Data` class** — a wrapper that supports arithmetic on nested numeric
  structures (scalars, lists, NumPy arrays).
- **Value display** — formatting of `(value, error)` pairs for publication.
- **Context managers** — temporary override of global jackknife and display
  settings.

```python
import qlat_utils as q
import numpy as np

avg, err = q.avg_err([1.0, 1.1, 0.9, 1.05])
print(q.show_val_err((avg, err)))
```

---

## Physical Constants

| Name | Value | Description |
|---|---|---|
| `alpha_qed` | `1 / 137.035999084` | Fine-structure constant |
| `fminv_gev` | `0.197326979` | Conversion factor: hbar*c / (1 fm * 1 GeV) |

```python
import qlat_utils as q
print(q.alpha_qed)   # 0.0072973525693...
print(q.fminv_gev)   # 0.197326979
```

---

## Type Tuples

Module-level tuples used for type-checking throughout the library. Extended
types (`float128`, `complex256`) are included when the platform supports them.

| Name | Contents |
|---|---|
| `float_types` | `float`, `np.float32`, `np.float64` (plus `np.float128` if available) |
| `complex_types` | `complex`, `np.complex64`, `np.complex128` (plus `np.complex256` if available) |
| `int_types` | `int`, `np.int32`, `np.int64` |
| `real_types` | `float_types + int_types` |
| `number_types` | `real_types + complex_types` |

---

## Interpolation

### `interp_i_arr(data_x_arr, x_arr)`

Return index array `i_arr` such that `q.interp(data_x_arr, i_arr)` is
approximately `x_arr`. Useful for mapping x-coordinates to fractional indices.

| Parameter | Type | Description |
|---|---|---|
| `data_x_arr` | array-like | Known x-values (must be monotonic) |
| `x_arr` | float or array-like | Target x-values |

### `interp(data_arr, i_arr, axis=-1)`

Return approximately `data_arr[..., i_arr]` using linear interpolation.
The index `i_arr` may be non-integer (fractional indices are interpolated
between adjacent elements).

| Parameter | Type | Description |
|---|---|---|
| `data_arr` | array-like | Source data |
| `i_arr` | float or 1-D array | Fractional index or indices |
| `axis` | int | Axis along which to interpolate (default `-1`) |

### `interp_x(data_arr, data_x_arr, x_arr, axis=-1)`

Interpolate `data_arr` at arbitrary x-values. Combines `interp_i_arr` and
`interp`.

| Parameter | Type | Description |
|---|---|---|
| `data_arr` | array-like | Source data |
| `data_x_arr` | array-like | x-values for `data_arr`; shape must be `(data_arr.shape[axis],)` |
| `x_arr` | float or 1-D array | Target x-values |
| `axis` | int | Axis along which to interpolate (default `-1`) |

### `get_threshold_idx(arr, threshold)`

Return the fractional index `x` such that `interp(arr, [x])` is approximately
`threshold`. Uses binary search on a 1-D array.

### `get_threshold_i_arr(data_arr, threshold_arr, axis=-1)`

Broadcast version of `get_threshold_idx` over an array. Returns an index array
where each entry satisfies the threshold condition along the given axis.

### `get_threshold_x_arr(data_arr, data_x_arr, threshold_arr, axis=-1)`

Like `get_threshold_i_arr`, but returns x-values instead of indices.

---

## Data Wrapper Class

### `class Data`

A wrapper around numeric values that supports arithmetic operations on nested
structures (scalars, lists, NumPy arrays, `LatData`).

**Supported value types:** numeric scalars, `numpy.ndarray`, `q.LatData`, and
`list` (element-wise operations).

```python
import qlat_utils as q

d1 = q.Data([1.0, 2.0, 3.0])
d2 = q.Data([0.5, 0.5, 0.5])
d3 = d1 + d2       # Data([1.5, 2.5, 3.5])
d4 = d1 * 2.0      # Data([2.0, 4.0, 6.0])
d5 = -d1           # Data([-1.0, -2.0, -3.0])
```

| Method | Description |
|---|---|
| `get_val()` | Return the wrapped value |
| `qnorm()` | Return the squared norm |
| `glb_sum()` | MPI global sum (requires `qlat`) |
| `__add__`, `__radd__` | Addition |
| `__sub__`, `__rsub__` | Subtraction |
| `__mul__`, `__rmul__` | Scalar or element-wise multiplication |
| `__neg__`, `__pos__` | Unary negation and identity |
| `__copy__`, `__deepcopy__` | Copy support |

---

## Basic Statistics

### `check_zero(x)`

Return `True` if `x` is a real type and equals zero.

### `qnorm(x)`

Return the squared norm of `x`. For scalars: `x*x`. For complex:
`re^2 + im^2`. For arrays: `abs(vdot(x, x))`. For lists/tuples: sum of
`qnorm` of each element.

```python
q.qnorm(2)          # 4
q.qnorm(1 + 2j)     # 5  (1*1 + 2*2)
```

### `average(data_list)`

Return the arithmetic mean of `data_list`.

### `average_ignore_nan(value_arr_list)`

Return element-wise average across a list of NumPy arrays, ignoring `NaN`
values. Returns `NaN` for elements where all inputs are `NaN`.

### `block_data(data_list, block_size, is_overlapping=True)`

Return a list of block averages. If `is_overlapping` is `True` (default),
blocks overlap by `block_size - 1` entries.

### `avg_err(data_list, *, eps=1, block_size=1)`

Compute `(avg, err)` of `data_list` using blocking. The error estimate is:

$$\text{err} = |\text{eps}| \sqrt{\frac{\text{block\_size}}{N - \text{block\_size}}} \cdot \text{fsqrt}\big(\text{avg}\big[(d_i - \text{avg})^2\big]\big)$$

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data_list` | list | — | Data values |
| `eps` | float | `1` | Additional scaling factor for error |
| `block_size` | int | `1` | Blocking size |

Returns `(avg, err)` where both have the same type as the data.

### `partial_sum(x, *, is_half_last=False)`

Modify `x` in-place to its cumulative (partial) sum, preserving length. If
`is_half_last` is `True`, each entry becomes the average of the current and
previous partial sums (trapezoidal rule). Works for 1-D and 2-D arrays.

### `fsqr(data)` / `fsqrt(data)`

Component-wise square and square root. For complex types, real and imaginary
parts are processed separately: `fsqr(a + bi) = a^2 + b^2 i`, `fsqrt(a + bi) =
sqrt(a) + sqrt(b) i`. Supports scalars, `Data`, and NumPy arrays.

### `err_sum(*vs)`

Return the quadrature sum of errors: `sqrt(sum(fsqr(v_i)))`.

```python
q.err_sum(1.4, 2.1, 1.0)  # 2.7147743920996454
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
| `n_rand_sample` | `1024` | Number of random samples (RJK only) |
| `is_normalizing_rand_sample` | `False` | Normalize random vectors (RJK only) |
| `is_hash_jk_idx` | `True` | Hash unknown jk indices (super only) |
| `jk_idx_hash_size` | `1024` | Hash table size (super only) |
| `block_size` | `1` | Default blocking size |
| `block_size_dict` | `{}` | Per-`job_tag` blocking sizes |
| `rng_state` | `RngState("rejk")` | RNG state |
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
summed with `Reduce_scatter` and `avg` is added to the samples.

### `get_gathered_jk_arr(jk_local, comm, num_node)`

Return the complete result from the local parts `jk_local`, gathered with
`Allgatherv` in the order of the nodes.

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

## Value Display

### `show_val(val, *, is_latex=True, num_float_digit=None, num_exp_digit=None, exponent=None)`

Format a single numeric value for display.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `val` | int/float | — | Value to format |
| `is_latex` | bool/None | `True` | Use LaTeX exponent notation |
| `num_float_digit` | int/bool/None | `None` | Number of decimal digits (auto if `None`) |
| `num_exp_digit` | int/bool/None | `None` | Significant digits in scientific notation |
| `exponent` | int/None | `None` | Force a specific exponent |

```python
q.show_val(0.00123)                     # "1.23 \\times 10^{-3}"
q.show_val(0.00123, is_latex=False)     # "1.23E-3"
q.show_val(1234.0)                      # "1234.0"
```

### `show_val_err(val_err, *, is_latex=True, num_float_digit=None, num_exp_digit=None, exponent=None)`

Format a `(value, error)` pair. Error is shown in parentheses. If `val_err` is
a single number, it is formatted as a plain value.

```python
q.show_val_err((1.12e16, 12e6))                         # auto scientific notation
q.show_val_err((1.12e16, 12e7), exponent=10)             # force exponent
q.show_val_err((1.12e16, 12e7), exponent=10, is_latex=False)  # "1.12000(120)E10"
```

---

## Context Managers

### `class NewDictValues(dictionary, **kwargs)`

Context manager that temporarily overrides keys in `dictionary` and restores
them on exit.

### `class JkKwargs(**kwargs)`

Context manager that temporarily overrides `default_g_jk_kwargs`.

```python
with q.JkKwargs(n_rand_sample=2048, block_size=10):
    jk_arr = q.g_mk_jk(data_list, jk_idx_list)
```

### `class ShowKwargs(**kwargs)`

Context manager that temporarily overrides `default_show_val_kwargs`.

```python
with q.ShowKwargs(is_latex=False, exponent=-10):
    print(q.show_val_err((1.23e-10, 0.05e-10)))
```

---

## Examples

### Interpolation

```python
import qlat_utils as q
import numpy as np

# Interpolate data at fractional indices
data = np.array([10.0, 20.0, 30.0, 40.0])
result = q.interp(data, 1.5)          # 25.0 (midpoint between 20 and 30)
result_arr = q.interp(data, [0.5, 1.5, 2.5])  # [15.0, 25.0, 35.0]

# Interpolate with explicit x-coordinates
x_data = np.array([0.0, 1.0, 2.0, 3.0])
y_data = np.array([0.0, 1.0, 4.0, 9.0])
x_new = np.array([0.5, 1.5, 2.5])
y_new = q.interp_x(y_data, x_data, x_new)  # interpolated y-values
```

### Basic Error Estimation

```python
import qlat_utils as q
import numpy as np

# Generate correlated data
data = [1.0 + 0.1 * np.random.randn() for _ in range(100)]

# Simple average and error
avg, err = q.avg_err(data)
print(f"avg = {avg:.4f}, err = {err:.4f}")

# With blocking to reduce autocorrelation
avg_b, err_b = q.avg_err(data, block_size=5)
print(f"avg = {avg_b:.4f}, err = {err_b:.4f}")
```

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

### Formatting Values

```python
import qlat_utils as q

# Format a single value
print(q.show_val(0.00123))                          # "1.23 \times 10^{-3}"
print(q.show_val(0.00123, is_latex=False))          # "1.23E-3"

# Format value with error
print(q.show_val_err((1.12e16, 12e6)))              # auto-notation
print(q.show_val_err((1.12e16, 12e7), exponent=10)) # "1.1200(12) \times 10^{10}"
```

### Context Managers

```python
import qlat_utils as q

# Temporarily change jackknife settings
with q.JkKwargs(n_rand_sample=2048, block_size=10):
    jk_arr = q.g_mk_jk(data_list, jk_idx_list)
    avg, err = q.g_jk_avg_err(jk_arr)

# Temporarily change display settings
with q.ShowKwargs(is_latex=False):
    print(q.show_val_err((3.14, 0.01)))
```

### Data Wrapper

```python
import qlat_utils as q

d1 = q.Data([1.0, 2.0, 3.0])
d2 = q.Data([0.1, 0.2, 0.3])

d3 = d1 + d2        # Data([1.1, 2.2, 3.3])
d4 = d1 * 2.0       # Data([2.0, 4.0, 6.0])
d5 = d1 - d2        # Data([0.9, 1.8, 2.7])
norm = d1.qnorm()   # 14.0 (1 + 4 + 9)
```
