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

The formula is:

$$jk\_arr[i] = \text{avg} + \sum_{j=1}^{N} \frac{-\text{eps}}{\sqrt{N(N - b(i,j))}} r_{i,j} (d_j - \text{avg})$$

where $r_{i,j} \sim \mathcal{N}(0, 1)$ and $b(i,j)$ is the block size.

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

### `g_mk_jk(data_list, jk_idx_list, *, avg=None, ...)`

Create a (randomized) super-jackknife dataset from un-jackknifed data.
Dispatches to `sjackknife` or `rjackknife` based on `jk_type`.

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
