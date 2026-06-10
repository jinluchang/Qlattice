# `qlat_utils.q_hlt_reconstruction` — HLT Spectral Reconstruction

Source: `qlat-utils/qlat_utils/q_hlt_reconstruction.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Parameters](#parameters)
   - [mk_hlt_params](#mk_hlt_params)
   - [get_f_e_weight_log](#get_f_e_weight_log)
3. [Core Transform](#core-transform)
   - [delta_from_g](#delta_from_g)
4. [Summation-Based Functions](#summation-based-functions)
   - [aa_from_g_via_sum](#aa_from_g_via_sum)
   - [normalization_constraint_via_sum](#normalization_constraint_via_sum)
   - [ww_from_g_via_sum](#ww_from_g_via_sum)
   - [mk_g_t_arr_via_sum](#mk_g_t_arr_via_sum)
5. [Integration-Based Functions](#integration-based-functions)
   - [aa_from_g](#aa_from_g)
   - [normalization_constraint](#normalization_constraint)
   - [ww_from_g](#ww_from_g)
   - [mk_g_t_arr](#mk_g_t_arr)
6. [Helper Functions](#helper-functions)
7. [Examples](#examples)

---

## Overview

The `qlat_utils.q_hlt_reconstruction` module implements the Hypercubic Laplacian
Transform (HLT) method for spectral density reconstruction. Given a target spectral
function `f_delta_target(e)` and noisy correlator data at discrete time slices, the
module finds a set of filter coefficients `g_t_arr` that reconstruct the spectral
density while minimizing a cost function that balances fidelity to the target and
regularization via a covariance penalty.

Two backends are provided:

- **Summation-based** (`_via_sum` suffix): evaluates integrals as discrete sums over
  a provided energy array `e_arr`. Faster but approximate.
- **Integration-based** (no suffix): uses `scipy.integrate.quad` for continuous
  energy integration. More accurate but slower.

Both backends use **JAX** for automatic differentiation (gradient computation via
`jax.value_and_grad`) and JIT compilation, enabling efficient optimization of the
filter coefficients.

The typical workflow is:

1. Create a parameter dictionary with `mk_hlt_params()`.
2. Set `f_delta_target`, `t_arr`, `cov`, and optionally `e_arr`.
3. Call `mk_g_t_arr(params)` or `mk_g_t_arr_via_sum(params)` to obtain optimal
   filter coefficients.
4. Use the coefficients to reconstruct the spectral density.

---

## Parameters

### `mk_hlt_params()`

Create and return a default HLT parameter dictionary.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `f_delta_target` | `callable` | `None` | Target spectral function `f(e)`. Required. |
| `t_arr` | `ndarray` | `None` | Array of discrete time values. Required. |
| `cov` | `float \| ndarray` | `None` | Covariance: scalar, 1D `(t_size,)`, or 2D `(t_size, t_size)`. Required. |
| `e_arr` | `ndarray` | `None` | Energy grid for summation-based functions. Required for `_via_sum` variants. |
| `e0` | `float` | `0.0` | Lower bound of energy integration. |
| `ee_max` | `float` | `np.inf` | Upper bound of energy integration (integration-based only). |
| `lambda` | `float` | `1.0` | Regularization weight for the covariance penalty term. |
| `alpha` | `float` | `-0.01` | Default exponent coefficient for the energy weight function. Overridden by `f_e_weight_log`. |
| `f_e_weight_log` | `callable \| None` | `None` | Custom log-weight function `f(e)`. If `None`, uses `alpha * e`. |
| `tt_size` | `float \| None` | `None` | Temporal extent for all-around-the-world (ATW) contributions. `None` disables ATW. |
| `atw_factor` | `float` | `1.0` | Multiplicative factor for ATW terms. Only effective when `tt_size` is not `None`. |
| `minimization_iter_max` | `int` | `1` | Number of minimization iterations. |
| `g_t_arr_init` | `ndarray \| None` | `None` | Initial guess for filter coefficients. If `None`, zeros are used. |
| `does_have_constraint` | `bool` | `False` | Whether to enforce a normalization constraint on the spectral density. |

### `get_f_e_weight_log(params)`

Return the log-weight function for the energy weight `exp(f(e))`.

If `params["f_e_weight_log"]` is set, returns that callable. Otherwise returns a
linear function `f(e) = alpha * e` using the `alpha` parameter.

---

## Core Transform

### `delta_from_g(g_t_arr, t_arr, e_arr)`

Compute the reconstructed spectral delta from filter coefficients.

**Parameters:**
- `g_t_arr` — Filter coefficients, shape `(..., t_size,)`.
- `t_arr` — Time array, shape `(t_size,)`.
- `e_arr` — Energy array, shape `(n_energies,)`.

**Returns:** Spectral delta, shape `(..., n_energies,)`.

Computes: `delta(e) = sum_t exp(-e * t) * g_t`

This is JIT-compiled with `@jax.jit`.

---

## Summation-Based Functions

These functions approximate integrals as discrete sums over the provided `e_arr`
energy grid. They are faster than the integration-based variants and suitable for
iterative optimization.

### `aa_from_g_via_sum(g_t_arr, params)`

Compute the chi-squared-like fidelity term using discrete summation.

**Parameters:**
- `g_t_arr` — Filter coefficients, shape `(t_size,)`.
- `params` — HLT parameter dictionary.

**Returns:** Scalar cost measuring the weighted squared difference between the
reconstructed and target spectral functions.

Evaluates: `sum_e w(e) * (delta(e) - delta_target(e))^2` where `w(e)` includes
trapezoidal quadrature weights and the energy weight function. Supports ATW
contributions when `tt_size` is set.

### `normalization_constraint_via_sum(g_t_arr, params)`

Enforce a normalization constraint via projection.

**Parameters:**
- `g_t_arr` — Filter coefficients, shape `(t_size,)`.
- `params` — HLT parameter dictionary.

**Returns:** Tuple `(new_g_t_arr, constraint_penalty)`.

If `does_have_constraint` is `False`, returns the input unchanged with zero
penalty. Otherwise, projects `g_t_arr` onto the subspace where the integrated
spectral density matches the target normalization, returning the projected
coefficients and a penalty proportional to the squared deviation.

### `ww_from_g_via_sum(g_t_arr, params)`

Compute the total cost function (fidelity + regularization + constraint).

**Parameters:**
- `g_t_arr` — Filter coefficients, shape `(t_size,)`.
- `params` — HLT parameter dictionary.

**Returns:** Scalar total cost.

The cost is: `aa / aa_zero + lambda * cov_term + constraint_penalty`, where
`aa_zero` is the fidelity of a zero filter (normalization baseline).

### `ww_from_g_wgrad_via_sum`

`jax.value_and_grad` of `ww_from_g_via_sum`. Returns `(value, gradient)`.

### `mk_g_t_arr_via_sum(params)`

Run the optimization to find optimal filter coefficients using summation.

**Parameters:**
- `params` — HLT parameter dictionary.

**Returns:** Optimized `g_t_arr` array.

Uses `q.q_fit_corr.minimize_scipy` for `minimization_iter_max` iterations,
starting from `g_t_arr_init` (or zeros). Decorated with `@q.timer`.

---

## Integration-Based Functions

These functions use `scipy.integrate.quad` for continuous energy integration. They
are more accurate but computationally heavier, making them better suited for final
results or small problems.

### `aa_from_g(g_t_arr, params)`

Compute the fidelity term using continuous integration.

**Parameters:**
- `g_t_arr` — Filter coefficients, shape `(t_size,)`.
- `params` — HLT parameter dictionary.

**Returns:** Scalar cost.

Internally builds and caches `aa_mat`, `f_vec`, and `aa_const` in `params` for
efficiency. The cost is: `aa_const + f_vec . g + g^T aa_mat g`.

### `normalization_constraint(g_t_arr, params)`

Enforce a normalization constraint using continuous integration.

**Parameters:**
- `g_t_arr` — Filter coefficients, shape `(t_size,)`.
- `params` — HLT parameter dictionary.

**Returns:** Tuple `(new_g_t_arr, constraint_penalty)`.

Same logic as `normalization_constraint_via_sum` but uses `build_hlt_fc_vec` for
continuous integration of the constraint vector.

### `ww_from_g(g_t_arr, params)`

Compute the total cost function using continuous integration.

**Parameters:**
- `g_t_arr` — Filter coefficients, shape `(t_size,)`.
- `params` — HLT parameter dictionary.

**Returns:** Scalar total cost.

Same structure as `ww_from_g_via_sum`: `aa / aa_zero + lambda * cov_term +
constraint_penalty`.

### `ww_from_g_wgrad`

`jax.value_and_grad` of `ww_from_g`. Returns `(value, gradient)`.

### `mk_g_t_arr(params)`

Run the optimization to find optimal filter coefficients using continuous
integration.

**Parameters:**
- `params` — HLT parameter dictionary.

**Returns:** Optimized `g_t_arr` array.

Same iterative procedure as `mk_g_t_arr_via_sum` but uses the integration-based
cost and gradient. Decorated with `@q.timer`.

---

## Helper Functions

### `get_cov_term(g_t_arr, cov)`

Compute the covariance regularization term `g^T Cov g`.

**Parameters:**
- `g_t_arr` — Filter coefficients, shape `(t_size,)`.
- `cov` — Covariance: 1D `(t_size,)` for diagonal or 2D `(t_size, t_size)` for full.

**Returns:** Scalar `sum(cov * g * g)` (diagonal) or `sum(cov * g_i * g_j)` (full).

### `build_hlt_aa_mat(params)`

Build and cache the quadratic coefficient matrix for the integration-based cost.
Matrix elements: `A[i,j] = integral exp(w(e) - e*(t_i + t_j)) de` with ATW terms.
Result is stored in `params["aa_mat"]`.

### `build_hlt_f_vec(params)`

Build and cache the linear coefficient vector for the integration-based cost.
`f[i] = integral -2 * delta_target(e) * exp(w(e) - e * t_i) de` with ATW terms.
Result is stored in `params["f_vec"]`.

### `build_hlt_fc_vec(params)`

Build and cache the normalization constraint vector. Returns `(norm, fc_vec)` where
`norm = integral delta_target(e) de` and the constraint is `norm == fc_vec . g`.
Result is stored in `params["fc_vec"]`.

### `build_hlt_aa_const(params)`

Build and cache the constant term in the integration-based cost.
`C = integral exp(w(e)) * delta_target(e)^2 de`.
Result is stored in `params["aa_const"]`.

---

## Examples

### Basic HLT Setup

```python
import qlat_utils as q
import qlat_utils.q_hlt_reconstruction
import numpy as np

# Create default parameters
params = q.q_hlt_reconstruction.mk_hlt_params()

# Set up a simple problem
t_size = 8
t_arr = np.arange(1, t_size, dtype=np.float64)
params["t_arr"] = t_arr
params["cov"] = 0.01 * np.ones(len(t_arr))
params["f_delta_target"] = lambda e: np.exp(-2.0 * e)

# Run optimization (summation-based)
e_arr = np.linspace(0.01, 5.0, 100)
params["e_arr"] = e_arr
g_t_arr = q.q_hlt_reconstruction.mk_g_t_arr_via_sum(params)
print(f"Filter coefficients: {g_t_arr}")
```

### Using Integration-Based Reconstruction

```python
import qlat_utils as q
import qlat_utils.q_hlt_reconstruction
import numpy as np

params = q.q_hlt_reconstruction.mk_hlt_params()
params["t_arr"] = np.arange(1, 8, dtype=np.float64)
params["cov"] = 0.01 * np.ones(7)
params["f_delta_target"] = lambda e: np.exp(-2.0 * e)
params["ee_max"] = 10.0

g_t_arr = q.q_hlt_reconstruction.mk_g_t_arr(params)
print(f"Filter coefficients: {g_t_arr}")
```

### With All-Around-the-World Contributions

```python
import qlat_utils as q
import qlat_utils.q_hlt_reconstruction
import numpy as np

params = q.q_hlt_reconstruction.mk_hlt_params()
params["t_arr"] = np.arange(1, 16, dtype=np.float64)
params["cov"] = 0.01 * np.ones(15)
params["f_delta_target"] = lambda e: np.exp(-3.0 * e)
params["tt_size"] = 32.0
params["atw_factor"] = 1.0
params["does_have_constraint"] = True

e_arr = np.linspace(0.01, 5.0, 200)
params["e_arr"] = e_arr
g_t_arr = q.q_hlt_reconstruction.mk_g_t_arr_via_sum(params)
print(f"Filter coefficients: {g_t_arr}")
```

### Evaluating the Cost Function

```python
import qlat_utils as q
import qlat_utils.q_hlt_reconstruction
import numpy as np

params = q.q_hlt_reconstruction.mk_hlt_params()
params["t_arr"] = np.arange(1, 8, dtype=np.float64)
params["cov"] = 0.01 * np.ones(7)
params["f_delta_target"] = lambda e: np.exp(-2.0 * e)

e_arr = np.linspace(0.01, 5.0, 100)
params["e_arr"] = e_arr

# Evaluate cost and gradient at zero filter
g_t_arr = np.zeros(7)
cost, grad = q.q_hlt_reconstruction.ww_from_g_wgrad_via_sum(g_t_arr, params)
print(f"Cost: {cost}")
print(f"Gradient: {grad}")
```
