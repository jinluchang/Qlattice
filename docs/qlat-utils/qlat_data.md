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
7. [Value Display](#value-display)
8. [Context Managers](#context-managers)
9. [Examples](#examples)

---

## Overview

The `qlat_utils.data` module provides the data analysis tools that are shared by
the other modules. It includes:

- **Interpolation** — linear interpolation for arrays with fractional indices.
- **Basic statistics** — averaging, blocking and the simple `avg_err` error
  estimation.
- **The `Data` class** — a wrapper that supports arithmetic on nested numeric
  structures (scalars, lists, NumPy arrays).
- **Value display** — formatting of `(value, error)` pairs for publication.
- **Generic helpers** — `use_kwargs`, the numeric type tuples and
  `NewDictValues`.

Jackknife resampling lives in `qlat_utils.jackknife_utils` (`q.jackknife`,
`q.g_mk_jk`, `q.g_jk_avg_err`, ...).

```python
import qlat_utils as q

avg, err = q.avg_err([1.0, 1.1, 0.9, 1.05])
print(q.show_val_err((avg, err)))
```

---

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

# Temporarily change the display settings.
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
