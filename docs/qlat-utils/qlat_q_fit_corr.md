# `qlat_utils.q_fit_corr` — Correlation-Matrix Fitting with Energy/Amplitude Parameterisation

Source: `qlat-utils/qlat_utils/q_fit_corr.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Correlation Model](#correlation-model)
3. [Data Generation](#data-generation)
4. [Objective Function Construction](#objective-function-construction)
5. [Minimisation Routines](#minimisation-routines)
6. [Multiprocessing Pool Helpers](#multiprocessing-pool-helpers)
7. [Jackknife Fitting: `fit_energy_amplitude`](#jackknife-fitting-fit_energy_amplitude)
8. [HMC Sampling](#hmc-sampling)
9. [Examples](#examples)

---

## Overview

`q_fit_corr` fits lattice-QCD two-point correlation matrices to a sum of
exponential decays parameterised by energies and operator amplitudes.
It uses **JAX** for automatic differentiation of the chi-squared objective and
**scipy** (or a built-in gradient-descent minimiser) for optimisation.
An optional Hamiltonian Monte Carlo (HMC) sampler is included for Bayesian
exploration of the posterior.

Key capabilities:

- Build synthetic correlation data with known parameters for testing.
- Construct JIT-compiled chi-squared functions with JAX.
- Run gradient-descent or L-BFGS-B minimisation with optional fixed/free
  parameter masks.
- Fit all jackknife samples in parallel via `multiprocessing`.
- Sample the posterior with a leapfrog HMC integrator with adaptive mass
  tuning.

---

## Correlation Model

The model correlation matrix is:

```
C_{i,j}(t) = sum_e  c[e,i] * c[e,j] * exp(-E[e] * (t + t_start[e]))
```

where `E[e]` are energies, `c[e,i]` are amplitudes, and `t_start[e]` are
optional per-state time offsets.

Optionally, around-the-world (ATW) periodicity effects are added:

```
C_{i,j}(t) += atw_factor[i] * atw_factor[j]
              * sum_e c[e,i] * c[e,j] * exp(-E[e] * (T-1-t + atw_t_start[e]))
```

Parameters are packed as `param_arr = [E_0, ..., E_{n-1}, c_00, c_01, ...]`.

---

## Data Generation

### `mk_data_set`

```python
q.mk_data_set(
    *, n_jk=10, n_ops=4, n_energies=4, t_size=4, sigma=0.1, rng=None
)
```

Generate synthetic jackknife correlation data with known ground-truth
parameters.

**Returns:** `(param_arr, jk_corr_data, corr_data_sigma)`

| Array | Shape | Description |
|---|---|---|
| `param_arr` | `(n_energies + n_energies*n_ops,)` | Packed energies and amplitudes |
| `jk_corr_data` | `(n_jk, n_ops, n_ops, t_size)` | Noisy jackknife samples |
| `corr_data_sigma` | `(n_ops, n_ops, t_size)` | Uniform noise sigma |

---

## Objective Function Construction

### `build_corr_from_param_arr`

```python
q.build_corr_from_param_arr(
    param_arr, *, n_ops, t_arr, t_start_arr=None,
    is_atw=False, atw_t_start_arr=None, atw_factor_arr=None,
)
```

Evaluate the model correlation matrix from packed parameters.

**Returns:** `corr` with shape `(n_ops, n_ops, len(t_arr))`.

### `mk_fcn`

```python
q.mk_fcn(
    corr_data, corr_data_sigma, t_start_arr, *,
    is_atw=False, atw_t_start_arr=None, atw_factor_arr=None,
    energy_minimum_arr=None, free_energy_idx_arr=None,
)
```

Build a JIT-compiled chi-squared function using JAX.

**Returns:** `fcn(param_arr, requires_grad=True)` that returns
`(chisq, grad)` or just `chisq`.

---

## Minimisation Routines

### `minimize`

```python
q.minimize(fcn, n_step=10, step_size=1e-2, *, param_arr)
```

Fixed-step-size gradient descent with early stopping when chi-squared
increases.

**Returns:** `(param_arr, n_steps_taken)`

### `adaptive_minimize`

```python
q.adaptive_minimize(
    fcn, step_size_list, n_step=10, max_total_steps=10000, *, param_arr,
)
```

Adaptive step-size gradient descent.  Cycles through `step_size_list`,
shrinking the step when the minimiser stalls and growing it when progress
is steady.

### `minimize_scipy`

```python
q.minimize_scipy(fcn, *, param_arr, fixed_param_mask=None, minimize_kwargs=None)
```

Scipy-based minimisation (default: L-BFGS-B) with optional fixed-parameter
mask.  Falls back to the initial parameters if the final objective is worse.

---

## Parameter Utilities

### `sort_param_arr_free_energy`

```python
q.sort_param_arr_free_energy(param_arr, n_ops, free_energy_idx_arr)
```

Sort free-energy states by ascending energy while keeping fixed states in
place.

### `apply_energy_minimum`

```python
q.apply_energy_minimum(param_arr, energy_minimum_arr=None, free_energy_idx_arr=None)
```

Clamp free energies to be at least `energy_minimum_arr` via the reflection
`E_min + |E - E_min|`.

---

## Multiprocessing Pool Helpers

### `mk_mp_pool` / `close_mp_pool`

```python
mp_pool = q.mk_mp_pool(n_proc=8)
# ... use mp_pool.imap ...
q.close_mp_pool(mp_pool)
```

Create a `multiprocessing.Pool` with spawn context and JAX pre-warming.

### `get_mp_pool_global` / `close_mp_pool_global`

```python
mp_pool = q.get_mp_pool_global(n_proc=8)
# ... use mp_pool.imap ...
q.close_mp_pool_global()
```

Module-level singleton pool.  Automatically re-created if `n_proc` changes.

---

## Jackknife Fitting: `fit_energy_amplitude`

```python
q.fit_energy_amplitude(
    jk_corr_data, *,
    t_start_data=0, t_start_fit=4, t_stop_fit=None,
    t_start_param=0, t_start_fcn=0,
    is_atw=False, atw_t_start_fcn=None, atw_factor=None,
    energy_minimum_arr=None,
    e_arr=None, c_arr=None,
    free_energy_idx_arr=None, fixed_coef_energy_idx_arr=None,
    n_step_mini_avg=10, n_step_mini_jk=5,
    minimize_kwargs=None, r_amp=1e-6,
    diag_err_scale_factor=1.0, off_diag_err_scale_factor=1.0,
    rng_seed_list=None, mp_pool=None,
)
```

Fit all jackknife samples of a correlation matrix.  The procedure is:

1. Normalise the correlation data by diagonal elements at `t=0`.
2. Run multiple random-start minimisations on the average data to find the
   best-fit parameters.
3. For each jackknife sample, minimise starting from the average-fit
   parameters (with optional random perturbation).

**Returns:** `res` dict with keys:

| Key | Shape | Description |
|---|---|---|
| `jk_chisq` | `(n_jk,)` | Chi-squared per jackknife sample |
| `jk_chisq_grad` | `(n_jk, n_params)` | Gradient of chi-squared |
| `jk_param_arr_for_scaled_corr` | `(n_jk, n_params)` | Parameters for scaled correlation |
| `jk_param_arr` | `(n_jk, n_params)` | Parameters in original normalisation |

---

## HMC Sampling

### `HmcParams`

Container for HMC trajectory state.  Key attributes:

| Attribute | Description |
|---|---|
| `traj` | Current trajectory number |
| `tau` | Molecular dynamics time |
| `n_step` | Number of leapfrog steps |
| `param_arr` | Current parameter vector |
| `hmc_mass_arr` | Per-parameter mass (inverse step-size scale) |
| `hmc_mass_adaptive_rate` | Rate for adaptive mass tuning |
| `temperature` | Heat-bath temperature |
| `delta_hh_history` | List of delta-H values |

### `hmc_traj`

```python
q.hmc_traj(fcn, hmc_params)
```

Run one HMC trajectory with leapfrog integration and Metropolis accept/reject.
Modifies `hmc_params` in place.

---

## Examples

### Fit synthetic data

```python
import qlat_utils as q
import numpy as np

param_true, jk_data, sigma = q.mk_data_set(
    n_jk=50, n_ops=3, n_energies=2, t_size=8, sigma=0.05,
)
e_arr = np.array([1.0, 2.0])
res = q.fit_energy_amplitude(
    jk_data, e_arr=e_arr,
    free_energy_idx_arr=np.array([0, 1]),
    n_step_mini_avg=5, n_step_mini_jk=3,
)
print("Fitted energies:", res["jk_param_arr"][:, :2].mean(axis=0))
```

### Custom minimisation with JAX

```python
import qlat_utils as q
import numpy as np

param_true, jk_data, sigma = q.mk_data_set(n_jk=5, t_size=6)
avg = jk_data.mean(axis=0)
fcn = q.mk_fcn(avg, sigma, np.zeros(4))
p0 = np.ones(4 * 5) * 0.5
p_opt, n_steps = q.minimize(fcn, n_step=50, step_size=0.01, param_arr=p0)
```
