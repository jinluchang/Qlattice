# `qlat_utils.utils` — General-Purpose Utility Functions

Source: `qlat-utils/qlat_utils/utils.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Environment & CLI Parsing](#environment--cli-parsing)
3. [Modular Arithmetic Helpers](#modular-arithmetic-helpers)
4. [Coordinate & Spatial Helpers](#coordinate--spatial-helpers)
5. [Lattice Physics Quantities](#lattice-physics-quantities)
6. [Interpolation](#interpolation)
7. [Data Signature & Test Verification](#data-signature--test-verification)
8. [Miscellaneous Utilities](#miscellaneous-utilities)
9. [Examples](#examples)

---

## Overview

The `qlat_utils.utils` module collects general-purpose functions used across
Qlattice packages. It re-exports symbols from `qlat_utils.timer`,
`qlat_utils.cache`, `qlat_utils.c`, and `qlat_utils.json`, so importing
`qlat_utils.utils` alone is sufficient for most scripts.

Key areas:

- **CLI argument and environment access** (`getenv`, `get_arg`, `get_option`, …)
- **Modular arithmetic** for lattice periodic boundary conditions (`rel_mod`, `rel_mod_arr`, …)
- **Spatial coordinate helpers** (`get_r_sq`, `get_r_limit`, `mk_r_list`, …)
- **Lattice physics quantities** (`epsilon_tensor`, `phat_sqr`, …)
- **Interpolation utilities** (`mk_interp_tuple`, `mk_r_sq_interp_idx_coef_list`)
- **Test result verification** (`json_results_append`, `check_log_json`)
- **Miscellaneous** (`lazy_call`, `get_fname`, `sqr`, `import_file`, `show_memory_usage`, …)

---

## Environment & CLI Parsing

### `getenv(*names, default=None)`

Return the value of the first defined environment variable among `names`, or
`default` if none is set. Logs the result via `displayln_info`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `*names` | `str` | — | One or more environment variable names |
| `default` | any | `None` | Fallback value |

### `get_arg(option, default=None, *, argv=None, is_removing_from_argv=False)`

Return the argument following `option` in `argv`, or `default` if `option` is
not found.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `option` | `str` | — | The CLI option to search for (e.g. `"--mass"`) |
| `default` | any | `None` | Returned when option is absent |
| `argv` | `list[str]` | `sys.argv` | Argument list to search |
| `is_removing_from_argv` | `bool` | `False` | Remove the option and arg from `argv` if `True` |

### `get_arg_list(option, *, argv=None, is_removing_from_argv=False)`

Collect **all** arguments following each occurrence of `option`. May appear
multiple times in `argv`.

Returns a `list[str]`.

### `get_option(option, *, argv=None, is_removing_from_argv=False)`

Return `True` if `option` is present in `argv`, `False` otherwise. Optionally
removes it.

### `get_all_arg_list(option, default=None, *, argv=None, is_removing_from_argv=False)`

Return all arguments that follow the **first** occurrence of `option` as a
list. Returns `default` if `option` is not found.

### `is_test()`

Return `True` if `--test` was passed on the command line (checked once at
module import time via `is_test_state`).

---

## Modular Arithmetic Helpers

These functions handle the periodic boundary conditions typical in lattice QCD,
returning results in the symmetric range `(-size/2, size/2]`.

### `rel_mod(x, size)`

Return `x % size` or `x % size - size` so the result is in `[0, size)` mapped
to `(-size/2, size/2]`.

### `rel_mod_sym(x, size)`

Like `rel_mod`, but returns `0` when `2 * x == size` (the exact midpoint).

### `rel_mod_arr(x, size)`

NumPy array version of `rel_mod`. `x` and `size` must be arrays of the same
shape.

### `rel_mod_sym_arr(x, size)`

NumPy array version of `rel_mod_sym`.

### `c_rel_mod(x, size)`

Component-wise `rel_mod` on sequences (lists/tuples) of the same length.

### `c_rel_mod_sqr(x, size)`

Sum of squared component-wise `rel_mod` values.

### `c_sqr(x)`

Sum of squares of all elements in sequence `x`.

---

## Coordinate & Spatial Helpers

### `parse_grid_coordinate_str(x_str)`

Parse a dot-separated coordinate string (e.g. `"2.3.1.0"`) into a
`Coordinate` object.

### `get_r_sq(x_rel)`

Return the spatial distance squared (as `int`) from the first three components
of `x_rel`.

### `get_r_limit(total_site)`

Return the maximum possible spatial `r` (as `float`) for a lattice with the
given `total_site` (a `Coordinate`). Computed as
`sqrt(sum((L_i / 2)^2))` over spatial dimensions.

### `mk_r_sq_list_3d(r_sq_limit)`

Return a sorted list of all unique `r^2 = x^2 + y^2 + z^2` values with
`r^2 <= r_sq_limit` for non-negative integers.

### `mk_r_sq_list(r_sq_limit, dimension="3D")`

Dispatch to `mk_r_sq_list_3d` for `"3D"` or return `range(0, r_sq_limit)` for
`"4D"` (Lagrange's four-square theorem).

### `mk_r_list(r_limit, *, r_all_limit=28.0, r_scaling_factor=5.0, dimension="3D")`

Generate a list of `r` values from `0` up to `r_limit`. Up to `r_all_limit`,
all discrete `r = sqrt(r_sq)` values are included. Beyond that, `r` values are
sampled at intervals of `1 / r_scaling_factor`.

---

## Lattice Physics Quantities

### `epsilon_tensor(i, j, k, l=3)`

Return the fully antisymmetric Levi-Civita tensor value for indices
`(i, j, k, l)`. `epsilon_tensor(0, 1, 2, 3) == 1`.

### `mk_epsilon_array()`

Construct and return the `4x4x4x4` NumPy `int8` array for the Levi-Civita
tensor. The module-level `epsilon_array` is precomputed from this.

### `phat_sqr(q, size)`

Compute the lattice momentum squared
`4 * sum_i sin^2(pi * (q_i % size_i) / size_i)`.

---

## Interpolation

### `mk_interp_tuple(x, x0, x1, x_idx)`

Return `(x_idx_low, x_idx_high, coef_low, coef_high)` for linear interpolation
of `x` between `x0` (at `x_idx`) and `x1` (at `x_idx + 1`).

### `mk_r_sq_interp_idx_coef_list(r_list)`

Build a list of interpolation tuples indexed by integer `r_sq`. Each entry
contains the two bounding `r_list` indices and linear interpolation
coefficients.

---

## Data Signature & Test Verification

### `get_data_sig_arr(x, rs, sig_len)`

Return a NumPy array of `sig_len` signature values extracted from the data `x`
(a `LatData`, `np.ndarray`, etc.) using `RngState` `rs`. The signature is
deterministic for a given data value and `rs`.

### `json_results_append(*args, json_results=None)`

Append a result tuple to `json_results` (defaults to `global_json_results`).
Used to collect test output for later comparison.

### `check_log_json(script_file, *, json_results=None, check_eps=1e-5)`

Compare the current `json_results` against a previously saved `.log.json` file.
Used by the CI system to detect regressions. Each result entry is
`(name, value, [check_eps])`.

### `global_json_results`

Module-level list used as the default accumulator for `json_results_append` and
`check_log_json`.

---

## Miscellaneous Utilities

### `lazy_call(f, *args, **kwargs)`

Return a zero-argument callable `get()` that evaluates `f(*args, **kwargs)` on
the first call and caches the result for subsequent calls.

### `get_fname()`

Return the function name of the caller (via `inspect.currentframe`).

### `sqr(x)`

Return `x * x`.

### `set_zero(x)`

Call `x.set_zero()`.

### `set_unit(x, coef=1.0)`

Call `x.set_unit(coef)`.

### `show(x)`

Call `x.show()`.

### `unitarize(x)`

Call `x.unitarize()`.

### `get_chunk_list(total_list, *, chunk_size=None, chunk_number=None, rng_state=None)`

Split `total_list` into chunks. Exactly one of `chunk_size` or
`chunk_number` must be provided. Optionally permute via `RngState`.

### `import_file(module_name, file_path)`

Import a Python module from an arbitrary file path and return it.

### `show_memory_usage()`

Print the current process RSS in GB (requires `psutil`).

### `displayln_info_malloc_stats()`

Print C++ malloc statistics on node 0.

---

## Examples

### Argument Parsing

```python
import qlat_utils as q

# Read an environment variable
mass = float(q.getenv("QLAT_MASS", "0.05"))

# Parse CLI: python script.py --mass 0.1 --verbose
mass = float(q.get_arg("--mass", "0.05"))
verbose = q.get_option("--verbose")
files = q.get_all_arg_list("--input")
```

### Modular Arithmetic

```python
import qlat_utils as q

# Periodic boundary with L=16
print(q.rel_mod(18, 16))       # -14
print(q.rel_mod_sym(8, 16))    # 0  (exact midpoint)

import numpy as np
x = np.array([18, 3, -1])
L = np.array([16, 16, 16])
print(q.rel_mod_arr(x, L))    # [-14   3  -1]
```

### Spatial Distance Lists

```python
import qlat_utils as q

r_list = q.mk_r_list(10.0)
print(r_list[:10])

r_sq_list = q.mk_r_sq_list(25, dimension="3D")
print(r_sq_list)
```

### Lattice Momentum

```python
import qlat_utils as q

size = [16, 16, 16, 32]
q_vec = [1, 0, 0, 0]
print(q.phat_sqr(q_vec, size))
```

### Lazy Evaluation

```python
import qlat_utils as q

expensive = q.lazy_call(int, "42")
print(expensive())  # evaluates on first call
print(expensive())  # returns cached result
```
