# `qlat_utils.q_fit_corr_2` — Correlation-Matrix Fitting with Eigenvalue/Coefficient Parameterisation

Source: `qlat-utils/qlat_utils/q_fit_corr_2.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Correlation Model](#correlation-model)
3. [Data Generation](#data-generation)
4. [Objective Function Construction](#objective-function-construction)
5. [Parameter Utilities](#parameter-utilities)
6. [Jackknife Fitting: `fit_eig_coef`](#jackknife-fitting-fit_eig_coef)
7. [Examples](#examples)

---

## Overview

`q_fit_corr_2` is a variant of `q_fit_corr` that parameterises the
correlation matrix using **eigenvalues** (propagator poles) and
**coefficients** rather than energies and amplitudes.  The key difference is
that the model uses power-law propagation `e^t` instead of exponential decay
`exp(-E*t)`, making it suitable for transfer-matrix eigenvalue fits where the
eigenvalues are extracted directly (e.g. from variational-basis analyses).

This module reuses `minimize_scipy`, `mk_mp_pool`, and `close_mp_pool` from
`q_fit_corr`.

Key capabilities:

- Build synthetic correlation data with eigenvalue/coefficient model.
- JIT-compiled chi-squared via JAX with optional eig-maximum constraints.
- Support for extra state signs (for negative-parity or alternating states).
- Around-the-world (ATW) periodicity effects.
- Parallel jackknife fitting with multiple random seeds.

---

## Correlation Model

The model correlation matrix is:

```
C_{i,j}(t) = sum_e  c[e,i] * c[e,j] * es[e]^(t - t_start[e]) * state_sign[e]
```

where `es[e]` are eigenvalues, `c[e,i]` are coefficients, and
`state_sign[e]` encodes sign factors for alternating states:

```
state_sign[e] = extra_state_sign[e] * (es[e] / |es[e]|)^((t_start[e] + extra_state_sign_t_start[e]) % 2)
```

Optionally, around-the-world effects are added:

```
C_{i,j}(t) += atw_factor[i] * atw_factor[j]
              * sum_e c[e,i] * c[e,j] * es[e]^(T - t - t_start[e]) * state_sign[e]
```

Parameters are packed as `param_arr = [es_0, ..., es_{n-1}, c_00, c_01, ...]`.

---

## Data Generation

### `mk_data_set`

```python
q.q_fit_corr_2.mk_data_set(
    *, n_jk=10, n_ops=4, n_eigs=4, t_arr=None, t_start_arr=None,
    extra_state_sign_arr=None, extra_state_sign_t_start_arr=None,
    t_size=None, atw_factor_arr=None, sigma=0.1, rng=None,
)
```

Generate synthetic jackknife correlation data with known eigenvalue/coefficient
ground truth.

**Returns:** `(param_arr, jk_corr_data, corr_data_sigma, t_arr)`

---

## Objective Function Construction

### `build_corr_from_param_arr`

```python
q.q_fit_corr_2.build_corr_from_param_arr(
    param_arr, *, t_arr, n_ops, n_eigs=None, t_start_arr=None,
    t_size=None, atw_factor_arr=None,
    extra_state_sign_arr=None, extra_state_sign_t_start_arr=None,
)
```

Evaluate the model correlation matrix from packed parameters.
Also accepts a 2-D `jk_param_arr` to build all jackknife samples at once.

**Returns:** `corr` with shape `(n_ops, n_ops, len(t_arr))` or
`(n_jk, n_ops, n_ops, len(t_arr))`.

### `mk_fcn`

```python
q.q_fit_corr_2.mk_fcn(
    corr_data, corr_data_sigma, t_arr, t_start_arr, *,
    extra_state_sign_arr=None, extra_state_sign_t_start_arr=None,
    t_size=None, atw_factor_arr=None,
    eig_maximum_arr=None, free_eig_idx_arr=None,
)
```

Build a JIT-compiled chi-squared function using JAX.

**Returns:** `fcn(param_arr, requires_grad=True)` that returns
`(chisq, grad)` or just `chisq`.

---

## Parameter Utilities

### `sort_param_arr_free_eig`

```python
q.q_fit_corr_2.sort_param_arr_free_eig(param_arr, n_ops, free_eig_idx_arr)
```

Sort free eigenvalue states by descending absolute value while keeping fixed
states in place.

### `apply_eig_maximum`

```python
q.q_fit_corr_2.apply_eig_maximum(param_arr, eig_maximum_arr=None, free_eig_idx_arr=None)
```

Clamp free eigenvalues to have magnitude at most `eig_maximum_arr` via the
reflection `max^2 / e` when `|e| > max`.

---

## Jackknife Fitting: `fit_eig_coef`

```python
q.q_fit_corr_2.fit_eig_coef(
    jk_corr_data, *,
    t_arr, e_arr, t_start_arr=None,
    extra_state_sign_arr=None, extra_state_sign_t_start_arr=None,
    t_size=None, atw_factor_arr=None,
    eig_maximum_arr=None,
    c_arr=None, op_norm_fac_arr=None,
    free_eig_idx_arr=None, fixed_coef_eig_idx_arr=None,
    n_step_mini_avg=10, n_step_mini_jk=5,
    minimize_kwargs=None, r_amp=1e-3,
    diag_err_scale_factor=1.0, off_diag_err_scale_factor=1.0,
    rng_seed_list=None, mp_pool=None,
)
```

Fit all jackknife samples using the eigenvalue/coefficient model.  The
procedure is:

1. Normalise correlation data by diagonal elements at `t[0]`.
2. Run multiple random-start minimisations on the average data.
3. For each jackknife sample, minimise starting from the best average-fit
   parameters.

**Returns:** `res` dict with keys:

| Key | Shape | Description |
|---|---|---|
| `jk_chisq` | `(n_jk,)` | Chi-squared per jackknife sample |
| `jk_param_arr` | `(n_jk, n_params)` | Parameters in original normalisation |
| `jk_param_grad_arr` | `(n_jk, n_params)` | Gradient of chi-squared |
| `jk_param_for_scaled_corr_arr` | `(n_jk, n_params)` | Parameters for scaled correlation |
| `jk_param_grad_for_scaled_corr_arr` | `(n_jk, n_params)` | Gradient for scaled correlation |
| `jk_e_arr` | `(n_jk, n_eigs)` | Eigenvalues |
| `jk_c_arr` | `(n_jk, n_eigs, n_ops)` | Coefficients |
| `jk_e_grad_arr` | `(n_jk, n_eigs)` | Eigenvalue gradients |
| `jk_c_grad_arr` | `(n_jk, n_eigs, n_ops)` | Coefficient gradients |
| `options` | dict | Echo of fit options (`t_arr`, `n_ops`, etc.) |

---

## Examples

### Fit synthetic eigenvalue data

```python
import qlat_utils as q
import numpy as np

q2 = q.q_fit_corr_2
t_arr = np.arange(6)
param_true, jk_data, sigma, t_arr = q2.mk_data_set(
    n_jk=50, n_ops=3, n_eigs=2, t_arr=t_arr, sigma=0.05,
)
e_arr = np.array([0.8, 0.5])
res = q2.fit_eig_coef(
    jk_data, t_arr=t_arr, e_arr=e_arr,
    free_eig_idx_arr=np.array([0, 1]),
    n_step_mini_avg=5, n_step_mini_jk=3,
)
print("Fitted eigenvalues:", res["jk_e_arr"].mean(axis=0))
```

### Use eig-maximum constraint

```python
import qlat_utils as q
import numpy as np

q2 = q.q_fit_corr_2
t_arr = np.arange(4)
param_true, jk_data, sigma, t_arr = q2.mk_data_set(
    n_jk=20, n_ops=2, n_eigs=2, t_arr=t_arr, sigma=0.1,
)
e_arr = np.array([0.9, 0.3])
res = q2.fit_eig_coef(
    jk_data, t_arr=t_arr, e_arr=e_arr,
    free_eig_idx_arr=np.array([0, 1]),
    eig_maximum_arr=np.array([1.0, 1.0]),
    n_step_mini_avg=3, n_step_mini_jk=2,
)
```
