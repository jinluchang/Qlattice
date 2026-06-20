# `qlat_scripts.v1.gen_data` — Propagator Generation and Field Selection Utilities

Source: `qlat/qlat_scripts/v1/gen_data.py`

> **Note:** Update this document when updating the source file.

## Overview

This module provides the core data generation pipeline for RBC/UKQCD lattice QCD calculations. It handles:

- Wall-source and point-source quark propagator inversions
- Smeared-source propagators
- Random U(1) volume-source propagators
- Field/point selection with importance sampling (weight-based)
- Sub-sampling of existing selections
- HVP (hadronic vacuum polarization) computation
- AMA (all-mode-accumulation) error correction support

All `run_*` functions follow a **lazy-evaluation** pattern: they return a callable `get_*` (wrapped with `@q.lazy_call`) that loads or computes data on first invocation and caches the result. If the output already exists on disk, the function returns the loader immediately without performing any computation.

## Terminology

| Term | Meaning |
|------|---------|
| `job_tag` | Ensemble identifier, e.g. `"24D"`, `"test-4nt8"` |
| `traj` | Trajectory number (integer) |
| `inv_type` | Quark flavor index: `0` = light, `1` = strange, `2` = charm |
| `inv_acc` | Inversion accuracy: `0` = sloppy, `1` = medium, `2` = exact |
| `gf` | Gauge field (`q.GaugeField`) |
| `gt` | Gauge transformation (`q.GaugeTransform`) |
| `eig` | Eigenvectors for deflation (from `run_eig`) |
| `psel` | Point selection (`q.PointsSelection`) |
| `fsel` | Field selection (`q.FieldSelection`) |
| `psel_prob` | Point selection with probability weights (`q.SelectedPointsRealD`) |
| `fsel_prob` | Field selection with probability weights (`q.SelectedFieldRealD`) |
| `wi` | Wall-source index list: list of `(idx, tslice, inv_type, inv_acc)` tuples |
| `sfw` | Selected-field writer (`q.open_fields(..., "a")`) |
| `qar_sp` | QAR archive for point-selected props (`q.open_qar_info(...)`) |

## Output Data Layout

All outputs are saved under `{job_tag}/...` with trajectory subdirectory `traj-{traj}/`. The typical naming conventions:

- **Full propagators**: `{job_tag}/prop-wsrc-full-{flavor}/traj-{traj}/`
- **Sparse wall-source props**: `{job_tag}/prop-wsrc-{flavor}/traj-{traj}/`
- **Point-source props**: `{job_tag}/prop-psrc-{flavor}/traj-{traj}/`
- **Smeared-source props**: `{job_tag}/prop-smear-{flavor}/traj-{traj}/`
- **Random U(1) props**: `{job_tag}/prop-rand-u1-{flavor}/traj-{traj}/`
- **Random U(1) sparse props**: `{job_tag}/prop-rand-u1-{type}-sparse-{flavor}/traj-{traj}/`
- **HVP fields**: `{job_tag}/hvp-psrc-{flavor}/traj-{traj}/`
- **Field selection weights**: `{job_tag}/field-selection-weight/traj-{traj}/`
- **Point-selected propagators**: `{job_tag}/psel-prop-{type}-{flavor}/traj-{traj}/`

where `{flavor}` is `"light"`, `"strange"`, or a name from `quark_flavor_list`.

Data is written atomically: files are first saved with `.acc` suffix, then renamed upon completion via `q.qrename_info`.

## Inverter Setup

### `run_get_inverter(job_tag, traj, *, inv_type, get_gf, get_gt=None, get_eig=None)`

Pre-computes and caches the inverter for all accuracy levels (`inv_acc` = 0, 1, 2) for a given quark flavor. Calls `ru.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)` internally.

**Parameters:**
- `get_gf` — callable returning gauge field (required)
- `get_gt` — callable returning gauge transform (optional, defaults to `None`)
- `get_eig` — callable returning eigenvectors (optional, defaults to `None`)

## Wall-Source Propagators

### `compute_prop_wsrc_1(job_tag, traj, *, gf, gt, eig, idx, tslice, inv_type, inv_acc)`

Core wall-source inversion. Creates a wall source at `tslice`, applies the inverter, and returns the solution propagator.

**Returns:** `q.Prop` — the solution field.

### `run_prop_wsrc_full(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt, get_wi)`

Computes **full** (un-sampled) wall-source propagators for all time slices and saves them along with their `qnorm_field` (used later for importance-sampled field selection). Skips if sparse wsrc data already exists.

**Output:** `{job_tag}/prop-wsrc-full-{flavor}/traj-{traj}/`

**Key parameters:**
- `get_wi` — callable returning the wall-source index list

### `run_prop_wsrc_sparse(job_tag, traj, *, inv_type, get_gf, get_gt, get_eig, get_psel, get_fsel, get_wi)`

Generates sparse (sampled) wall-source propagators by loading full wsrc data or performing on-the-fly inversions, then projecting onto `psel` and `fsel`. Saves both the selected-field propagator and the wall-sink projected propagator.

**Output:**
- `{job_tag}/prop-wsrc-{flavor}/traj-{traj}/` (selected-field data)
- `{job_tag}/psel-prop-wsrc-{flavor}/traj-{traj}/` (point-selected data + wall-sink)

**Behavior:** If full wsrc data is unavailable and `is_performing_inversion_if_no_full_prop_available` is `False` (default), prints a warning and returns.

## Field Selection Weight Computation

### `run_f_weight_from_wsrc_prop_full(job_tag, traj)`

Computes importance-sampling weights from full wall-source propagator norms. Returns `get_f_weight` callable that yields a `q.FieldRealD(geo, 1)` with per-site weights averaged around 1.

**Algorithm:**
1. Loads `qnorm_field` from full wsrc data for both light and strange quarks
2. Computes global-sum-per-tslice profiles
3. Averages weight profiles across time slices
4. Combines light (25%) and strange (25%) with uniform baseline (50%)

**Output:** `{job_tag}/field-selection-weight/traj-{traj}/weight.field`

**Returns:** `get_f_weight` callable, or `None` if data is unavailable.

### `run_f_weight_uniform(job_tag, traj)`

Alternative to `run_f_weight_from_wsrc_prop_full`: creates a uniform weight field (`f_weight.set_unit()`). Useful for testing or when full wsrc data is not yet available.

**Returns:** `get_f_weight` callable.

### `run_f_weight_load(job_tag, traj)`

Loads an existing weight field from disk. Raises `Exception` if the file does not exist.

**Returns:** `get_f_weight` callable.

## Random Field and Selection Probability

### `run_f_rand_01(job_tag, traj)`

Generates a reproducible uniform random field in [0, 1) used for stochastic selection. The RNG seed is derived from `get_job_seed(job_tag)` and the trajectory number.

**Output:** `{job_tag}/field-selection-weight/traj-{traj}/f-rand-01.field`

**Returns:** `get_f_rand_01` callable yielding `q.FieldRealD(geo, 1)`.

### `run_fsel_prob(job_tag, traj, *, get_f_rand_01, get_f_weight)`

Creates the field selection (stochastic spatial sampling) with probability weights.

**Algorithm:**
1. Selects sites where `f_weight * fsel_rate >= f_rand_01`
2. Saves the `FieldSelection` and a `SelectedFieldRealD` of selection probabilities

**Parameters:**
- `get_f_rand_01` — set to `None` to load existing data
- `get_f_weight` — set to `None` to load existing data

**Output:**
- `{job_tag}/field-selection/traj-{traj}.field`
- `{job_tag}/field-selection-weight/traj-{traj}/fsel-prob.sfield`

**Returns:** `get_fsel_prob` callable yielding `q.SelectedFieldRealD(fsel, 1)`.

### `run_psel_prob(job_tag, traj, *, get_f_rand_01, get_f_weight, tag=None)`

Creates the point selection with probability weights. Same algorithm as `run_fsel_prob` but for points.

**Parameters:**
- `tag` — optional tag for named selections (e.g. `"small"`, `"median"`, `"large"`)

**Output:**
- `{job_tag}/points-selection/traj-{traj}.lati` (default, no tag)
- `{job_tag}/field-selection-weight/traj-{traj}/psel-prob.lat` (default, no tag)
- `{job_tag}/psel_{tag}/traj-{traj}/psel.lati` (with tag)
- `{job_tag}/psel_{tag}/traj-{traj}/psel-prob.lat` (with tag)

**Returns:** `get_psel_prob` callable yielding `q.SelectedPointsRealD(psel, 1)`.

### `run_fsel_from_fsel_prob(get_fsel_prob)` / `run_psel_from_psel_prob(get_psel_prob)`

Convenience wrappers that extract `fsel` / `psel` from the probability-weighted objects.

**Returns:** `lambda: get_fsel_prob().fsel` / `lambda: get_psel_prob().psel`, or `None` if input is `None`.

## Sub-Sampling

### `run_fsel_prob_sub_sampling(job_tag, traj, *, sub_sampling_rate, get_fsel_prob, get_f_rand_01, get_f_weight)`

Creates a sub-sample of an existing field selection. Approximately `sub_sampling_rate` fraction of the original selection is kept.

**Parameters:**
- `sub_sampling_rate` — fraction in [0, 1]; `1.0` means complete sub-sampling
- `get_f_weight` — if `None`, uses `fsel_prob * sub_sampling_rate` as probability (not exactly equivalent to using `f_weight`)

**Returns:** `get_fsel_prob_sub` callable.

### `run_psel_prob_sub_sampling(job_tag, traj, *, sub_sampling_rate, get_psel_prob, get_f_rand_01, get_f_weight)`

Point-selection analogue of `run_fsel_prob_sub_sampling`.

If `get_param(job_tag, "use_simple_psel_sub_sampling", default=False)` is `True`, simply takes the first `sub_sampling_rate * original_num` points.

**Returns:** `get_psel_prob_sub` callable.

## Selection Splitting

### `run_psel_split(job_tag, traj, *, get_psel, num_piece)` / `run_fsel_split(job_tag, traj, *, get_fsel, num_piece)`

Splits a point/field selection into `num_piece` sub-selections with increased spatial separation (for independent measurement chunks). `num_piece` should be a power of 2.

**Output:** `{job_tag}/points-selection-split/traj-{traj}/num-piece-{num_piece}/` or `{job_tag}/field-selection-split/traj-{traj}/num-piece-{num_piece}/`

**Returns:** `get_psel_list` callable yielding `list[q.PointsSelection]` of length `num_piece`.

## Point-Source Propagators

### `compute_prop_psrc(job_tag, traj, xg_src, inv_type, inv_acc, *, idx, gf, gt, sfw, qar_sp, psel, fsel, f_rand_01, sfw_hvp, qar_hvp_ts, eig)`

Core point-source inversion at position `xg_src`. Performs AMA-style field selection for the solution: sites with large propagator norm (above `field_selection_fsel_psrc_prop_norm_threshold`) are probabilistically added to the field selection.

**Saves:**
- Selected-field propagator (in `sfw`)
- Point-selected propagator + wall-sink projection (in `qar_sp`)
- Additional `fsel-prob-psrc-prop` field for the adaptive selection
- HVP contraction (if `sfw_hvp` / `qar_hvp_ts` are provided)

### `run_prop_psrc(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt, get_psel, get_fsel, get_f_rand_01)`

Runs point-source propagator generation for all points in `psel` with AMA multi-accuracy. For each source point, sloppy inversion is always performed; medium and exact inversions are done stochastically based on `prob_acc_1_psrc` and `prob_acc_2_psrc`.

**Output:**
- `{job_tag}/prop-psrc-{flavor}/traj-{traj}/`
- `{job_tag}/psel-prop-psrc-{flavor}/traj-{traj}/`
- `{job_tag}/hvp-psrc-{flavor}/traj-{traj}/` (HVP fields)
- `{job_tag}/hvp-sum-tslice-psrc-{flavor}/traj-{traj}/` (HVP time-slice sums)

## HVP Computation

### `calc_hvp_sum_tslice(chvp_16)`

Computes the HVP summed over spatial slices in all 4 directions from a `chvp_16` field (produced by `q.contract_chvp_16`).

**Returns:** `ld_hvp_ts` — LatData with shape `[t_dir, t, mu, nu]` where `t_dir` ∈ {x, y, z, t}.

### `compute_hvp_average(job_tag, traj, *, inv_type, psel_prob, data_path, geo)`

Computes the AMA-corrected average HVP field from point-source HVP data, incorporating probability weights and source-position shifting.

**Returns:** `hvp_average` — `q.FieldComplexD(geo, 16)`.

### `run_hvp_average(job_tag, traj, *, inv_type, get_psel_prob)`

Wrapper that loads or computes the average HVP field.

**Output:** `{job_tag}/hvp-average/traj-{traj}/hvp_average_{flavor}.field`

**Returns:** `load` callable yielding `q.FieldComplexD(geo, 16)`.

## Random U(1) Volume-Source Propagators

### `run_field_rand_u1_dict(job_tag, traj)`

Generates reproducible random U(1) fields for both `fsel` and `psel` source types, plus their conjugates.

**Output:** `{job_tag}/field-rand-u1/traj-{traj}/`

**Returns:** `get_field_rand_u1_dict` callable yielding a dict with keys `"fsel-src"`, `"fsel-src-dag"`, `"psel-src"`, `"psel-src-dag"`, each mapping to a `q.FieldComplexD`.

### `run_prop_sparse_rand_u1_src(job_tag, traj, *, inv_type, get_gf, get_psel, get_fsel, get_field_rand_u1_dict, get_psel_list=None, get_fsel_psel_list=None, get_eig=None)`

Computes propagators from random U(1) sparse sources defined on `psel_list` or `fsel_psel_list` sub-selections. Supports both dagger and non-dagger inversions with AMA multi-accuracy.

**Parameters:**
- Set exactly one of `get_psel_list` or `get_fsel_psel_list`
- `get_psel_list` — callable returning `list[q.PointsSelection]` for psel-type sources
- `get_fsel_psel_list` — callable returning `list[q.PointsSelection]` for fsel-type sources

**Output:**
- `{job_tag}/prop-rand-u1-{type}-sparse-{flavor}/traj-{traj}/`
- `{job_tag}/psel-prop-rand-u1-{type}-sparse-{flavor}/traj-{traj}/`

### `run_prop_rand_u1(job_tag, traj, *, inv_type, get_gf, get_fsel, get_eig=None)`

Computes random U(1) volume-source propagators using `q.mk_rand_u1_prop`. Each source uses a unique RNG seed derived from `(job_seed, traj, idx_rand_u1)`.

**Parameters:**
- Number of sources controlled by `get_param(job_tag, "n_rand_u1_fsel")`

**Output:** `{job_tag}/prop-rand-u1-{flavor}/traj-{traj}/`

## Smeared-Source Propagators

### `run_prop_smear(job_tag, traj, *, inv_type, get_gf, get_gf_ape, get_eig, get_gt, get_psel, get_fsel, get_psel_smear, get_psel_smear_median)`

Generates propagators from APE-smeared sources. The smearing parameters are read from job parameters:
- `get_param(job_tag, "prop_smear_coef")`
- `get_param(job_tag, "prop_smear_step")`

Both the original and smeared-sink solutions are saved.

**Output:**
- `{job_tag}/prop-smear-{flavor}/traj-{traj}/`
- `{job_tag}/psel-prop-smear-{flavor}/traj-{traj}/`
- `{job_tag}/psel_smear_median-prop-smear-{flavor}/traj-{traj}/`

## Typical Workflow

A typical data generation script (see `examples-py-gpt/gpt-qlat-auto-simple.py`) follows this sequence:

```python
from qlat_scripts.v1 import *

# 1. Gauge field and transform
get_gf = run_gf(job_tag, traj_gf)
get_gt = run_gt(job_tag, traj_gf, get_gf)

# 2. Eigenvectors (for light quark deflation)
get_eig_light = run_eig(job_tag, traj_gf, get_gf)
get_eig_strange = run_eig_strange(job_tag, traj_gf, get_gf)

# 3. Wall-source index list
get_wi = run_wi(job_tag, traj)

# 4. Full wall-source propagators (for importance sampling weights)
run_prop_wsrc_full(job_tag, traj, inv_type=0, get_gf=get_gf,
                   get_eig=get_eig_light, get_gt=get_gt, get_wi=get_wi)
run_prop_wsrc_full(job_tag, traj, inv_type=1, get_gf=get_gf,
                   get_eig=get_eig_strange, get_gt=get_gt, get_wi=get_wi)

# 5. Field/point selection from wsrc propagator norms
get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj)
get_f_rand_01 = run_f_rand_01(job_tag, traj)
get_fsel_prob = run_fsel_prob(job_tag, traj,
                              get_f_rand_01=get_f_rand_01,
                              get_f_weight=get_f_weight)
get_psel_prob = run_psel_prob(job_tag, traj,
                              get_f_rand_01=get_f_rand_01,
                              get_f_weight=get_f_weight)
get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
get_psel = run_psel_from_psel_prob(get_psel_prob)

# 6. Sparse wall-source propagators (sampled from full data)
run_prop_wsrc_sparse(job_tag, traj, inv_type=0, get_gf=get_gf,
                     get_gt=get_gt, get_eig=get_eig_light,
                     get_psel=get_psel, get_fsel=get_fsel, get_wi=get_wi)

# 7. Point-source propagators (with AMA)
run_prop_psrc(job_tag, traj, inv_type=0, get_gf=get_gf,
              get_eig=get_eig_light, get_gt=get_gt,
              get_psel=get_psel, get_fsel=get_fsel,
              get_f_rand_01=get_f_rand_01)

# 8. Random U(1) propagators
run_prop_rand_u1(job_tag, traj, inv_type=0, get_gf=get_gf,
                 get_fsel=get_fsel, get_eig=get_eig_light)
```

### Alternative: Uniform Weights (Testing)

For testing or when full wsrc data is unavailable:

```python
get_f_weight = run_f_weight_uniform(job_tag, traj)
```

### Multi-Flavor with Extended Quark Masses

For ensembles with charm-like quarks (see `examples-py-gpt/gpt-qlat-data-gen-prop.py`):

```python
quark_mass_list = get_param(job_tag, "quark_mass_list")
for inv_type, mass in enumerate(quark_mass_list):
    if inv_type == 0:
        get_eig = get_eig_light
    else:
        get_eig = None
    run_prop_rand_vol_u1_src(job_tag, traj, inv_type=inv_type,
                             get_gf=get_gf, get_psel=get_psel,
                             get_fsel=get_fsel, get_eig=get_eig)
```

## Key Parameters

These parameters are read via `get_param(job_tag, ...)`:

| Parameter | Description |
|-----------|-------------|
| `total_site` | Lattice dimensions, e.g. `[24, 24, 24, 64]` |
| `fermion_params` | Fermion action parameters per `(inv_type, inv_acc)` |
| `cg_params-{inv_type}-{inv_acc}` | CG solver parameters (`maxiter`, `maxcycle`) |
| `field_selection_fsel_rate` | Fraction of sites in field selection |
| `field_selection_psel_rate` | Fraction of sites in point selection |
| `field_selection_fsel_psrc_prop_norm_threshold` | Threshold for adaptive psrc field selection |
| `prob_exact_wsrc` | Probability of exact-accuracy wall-source inversion |
| `prob_acc_1_psrc` | Probability of medium-accuracy point-source inversion |
| `prob_acc_2_psrc` | Probability of exact-accuracy point-source inversion |
| `prob_acc_1_rand_u1` | Probability of medium-accuracy rand U(1) inversion |
| `prob_acc_2_rand_u1` | Probability of exact-accuracy rand U(1) inversion |
| `n_rand_u1_fsel` | Number of random U(1) sources |
| `prop_smear_coef` | APE smearing coefficient |
| `prop_smear_step` | APE smearing steps |
| `quark_flavor_list` | List of flavor names, e.g. `["light", "strange", "charm-1"]` |
| `quark_mass_list` | List of quark masses |
| `lanc_params` | Lanczos parameters for eigenvector generation |
| `clanc_params` | Chebyshev-Lanczos parameters |

## Notes

- All `run_*` functions use `q.obtain_lock()` / `q.release_lock()` for file-based mutual exclusion in multi-node environments.
- AMA (all-mode-accumulation) is implemented stochastically: sloppy inversions are always performed, while medium and exact inversions are done with probabilities `prob_acc_1` and `prob_acc_2`. The probability is used later to reweight results.
- The `@q.lazy_call` decorator ensures each `get_*` callable executes its computation at most once.
- The `@q.timer` and `@q.timer_verbose` decorators provide automatic timing instrumentation.
- Functions decorated with `@q.timer(is_timer_fork=True)` support timer forking for parallel execution contexts.
