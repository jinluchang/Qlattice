# `qlat_scripts.v1.jobs` — Core Job Management

Source: `qlat/qlat_scripts/v1/jobs.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Path Management](#path-management)
3. [Job Checking](#job-checking)
4. [Parameter Management](#parameter-management)
5. [Gauge Fields](#gauge-fields)
6. [Wall-Source Info](#wall-source-info)
7. [Point Selection](#point-selection)
8. [Point-Source Info](#point-source-info)
9. [Point Distribution](#point-distribution)
10. [Field Selection](#field-selection)
11. [Smeared Selections](#smeared-selections)
12. [Gauge Smearing](#gauge-smearing)
13. [Eigenvectors](#eigenvectors)
14. [Radial Lists](#radial-lists)
15. [Examples](#examples)

---

## Overview

This module is the central job management layer for the `qlat_scripts.v1` pipeline. It provides:

- **Path management**: save/load path resolution with multiple search directories
- **Job scheduling**: `check_job` for determining if a trajectory needs processing
- **Gauge field loading**: `run_gf` for gauge configurations, `run_gt` for Coulomb gauge fixing
- **Source selection**: wall-source info, point-source info, point selection, field selection
- **Smeared selections**: APE-smeared point selections for extended operators
- **Gauge smearing**: APE and HYP smearing of gauge fields
- **Eigenvector computation**: Lanczos/chebyshev-Lanczos eigenvectors for deflation
- **Radial basis functions**: momentum-projected correlator support via `r_list`

All `run_*` functions follow a **lazy-evaluation** pattern: they return a callable `get_*` (wrapped with `@q.lazy_call`) that loads or computes data on first invocation and caches the result. If the output already exists on disk, the function returns the loader immediately.

---

## Path Management

### `save_path_default`

Default save directory: `"results"`.

### `load_path_list`

List of directories to search when loading data: `["results"]`.

### `get_save_path(fn)`

Returns `f"{save_path_default}/{fn}"`.

### `get_load_path(*fns)`

Searches `load_path_list` for the first existing file matching any of the given filenames. Supports nested lists/tuples. Uses `q.does_file_exist_qar_sync_node` for QAR-aware file existence checks.

---

## Job Checking

### `check_job(job_tag, traj, fns_produce, fns_need)`

Determines whether a trajectory needs processing.

**Returns:** `True` if the job should be run (outputs missing, inputs available). Returns `False` if all outputs exist (finished) or any required input is missing (unavailable).

Also calls `q.check_stop()` and `q.check_time_limit()` to support graceful termination.

---

## Parameter Management

### `run_params(job_tag)`

Saves ensemble parameters to `{job_tag}/params/version-{N}.json` and `.pickle`. Version numbers are zero-padded 10-digit integers. If the latest saved version matches the current parameters, no new version is written.

### `set_param(job_tag, tag, value)` / `get_param(job_tag, tag, default=None)`

Re-exported from `rbc_ukqcd_params`. Accesses the global `dict_params` dictionary.

### `get_job_seed(job_tag)`

Returns the RNG seed for the given ensemble (re-exported from `rbc_ukqcd_params`).

---

## Gauge Fields

### `mk_gf_fn_list(job_tag, traj)`

Returns a list of candidate filenames for gauge configurations, supporting multiple streaming offsets (up to `4 * 1000000`). Checks both native and IEEE64BIG formats across `configs/`, `configs-b/`, etc.

### `run_gf(job_tag, traj)`

Loads a gauge field. For `"test-*"` ensembles, generates a sample gauge field via `rup.mk_sample_gauge_field_v3`. Returns a lazy `get_gf` callable, or `None` if the configuration is not found.

### `run_gt(job_tag, traj, get_gf)`

Computes or loads the Coulomb gauge transformation. Uses `qlat_gpt.gauge_fix_coulomb`. Saves both CPS and native formats. Returns a lazy `get_gt` callable, or `None` if inputs are unavailable.

---

## Wall-Source Info

Wall-source info (`wi`) is a list of `[idx, tslice, inv_type, inv_acc]` tuples specifying which time slices get sloppy vs. exact inversions.

### `mk_rand_wall_src_info_n_exact(job_tag, traj, inv_type)`

Generates wall-source info with exactly `n_exact_wsrc` randomly chosen time slices at `inv_acc=2` (exact), all others at `inv_acc=1` (sloppy).

### `mk_rand_wall_src_info_prob(job_tag, traj, inv_type)`

Generates wall-source info where each time slice independently gets exact accuracy with probability `prob_exact_wsrc`.

### `mk_rand_wall_src_info(job_tag, traj, inv_type)`

Dispatches to `_n_exact` or `_prob` depending on whether `prob_exact_wsrc` is in the ensemble parameters.

### `get_prob_exact_wsrc(job_tag)`

Returns the probability of exact wall-source inversion. Computes from `n_exact_wsrc` if `prob_exact_wsrc` is not set.

### `save_wall_src_info(wi, path)` / `load_wall_src_info(path)`

Serialize/deserialize wall-source info to/from text files.

### `run_wi(job_tag, traj)`

Generates or loads wall-source info for both light and strange quarks. Returns a lazy callable yielding the combined list.

---

## Point Selection

### `get_n_points_psel(job_tag)`

Returns the number of points in the point selection. Uses `field_selection_psel_rate` if available, otherwise `n_points_psel`.

### `mk_rand_psel(job_tag, traj)`

Generates a random point selection with `n_points` points.

### `run_psel(job_tag, traj)`

Generates or loads a random point selection. Returns a lazy `get_psel` callable.

---

## Point-Source Info

Point-source info (`pi`) is a list of `[idx, xg, inv_type, inv_acc]` tuples specifying source positions and accuracy levels.

### `get_n_points_pi(job_tag, traj, inv_type, inv_acc)`

Returns the number of point sources for a given accuracy level.

### `mk_rand_point_src_info(job_tag, traj, psel)`

Generates point-source info grouping multiple accuracy levels per source position.

### `save_point_src_info(pi, path)` / `load_point_src_info(path)`

Serialize/deserialize point-source info.

### `run_pi(job_tag, traj, get_psel)`

Generates or loads point-source info. Returns a lazy callable.

---

## Point Distribution

### `load_point_distribution(job_tag)`

Loads relative coordinate probability distribution for point selections. Returns a dict mapping `(x, y, z, t)` tuples (with `x >= y >= z >= 0`) to probabilities.

### `classify_rel_coordinate(xg_rel_arrary, total_site_array)`

Classifies a relative coordinate into its canonical form with `x >= y >= z >= 0` and periodic boundary conditions.

### `get_point_xrel_prob(xg_rel_arrary, total_site_array, point_distribution, n_points)`

Returns the probability for a given relative coordinate. Falls back to uniform distribution if `point_distribution` is `None`.

---

## Field Selection

### `mk_rand_fsel(job_tag, traj, n_per_tslice)`

Generates a random field selection with `n_per_tslice` sites per time slice.

### `run_fsel(job_tag, traj)`

Generates or loads a field selection. The selection rate is controlled by `field_selection_fsel_rate`. Returns a lazy `get_fsel` callable.

### `mk_fselc(fsel, psel)`

Creates a combined field selection that is guaranteed to contain all points in `psel`.

### `run_fselc(job_tag, traj, get_fsel, get_psel)`

Lazy wrapper around `mk_fselc`.

---

## Smeared Selections

### `mk_rand_fsel_smear(job_tag, traj, n_per_tslice_smear)`

Generates a random selection for smeared-source operators.

### `run_psel_smear(job_tag, traj)`

Generates or loads a smeared point selection with uniform sampling per time slice. Returns a lazy callable.

### `run_psel_smear_median(job_tag, traj)`

Generates or loads a median-density smeared point selection. Returns a lazy callable.

---

## Gauge Smearing

### `run_gf_ape(job_tag, get_gf)`

Applies APE spatial smearing to the gauge field. Parameters `gf_ape_smear_coef` and `gf_ape_smear_step` from `dict_params`. Returns a lazy `get_gf_ape` callable.

### `run_gf_hyp(job_tag, get_gf)`

Applies HYP smearing to the gauge field. Parameters `gf_hyp_smear_step` from `get_param`. Uses fixed coefficients `(0.75, 0.6, 0.3)`. Returns a lazy `get_gf_hyp` callable.

---

## Eigenvectors

### `compute_eig(job_tag, gf, inv_type=0, inv_acc=0, *, path=None, pc_ne=None)`

Computes or loads eigenvectors via `ru.mk_ceig`. Saves with atomic rename (`.partial` suffix). Runs `test_eig` to validate. Returns a `get_eig` callable.

### `test_eig(job_tag, gf, eig, inv_type, *, pc_ne=None)`

Validates eigenvectors by comparing sloppy and medium accuracy inversions (with and without deflation) against a reference exact solution. Results are appended to `q.json_results_append`.

### `run_eig(job_tag, traj, get_gf, *, is_only_load=False)`

Loads or computes light-quark eigenvectors. Returns `None` if unavailable.

### `run_eig_strange(job_tag, traj, get_gf, *, is_only_load=False)`

Loads or computes strange-quark eigenvectors (chebyshev-Lanczos). Returns `lambda: None` if `clanc_params` are not configured.

---

## Radial Lists

### `get_r_list(job_tag)`

Returns the list of radial distances `r` for momentum-projected correlators. Cached via `@functools.lru_cache`.

### `get_r_sq_interp_idx_coef_list(job_tag)`

Returns interpolation coefficients `(r_idx_low, r_idx_high, coef_low, coef_high)` indexed by `r_sq`. Cached via `@functools.lru_cache`.

### `run_r_list(job_tag)`

Saves the radial list to `{job_tag}/r_list/r_list.lat`. Validates against existing data if present.

---

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])
import qlat_scripts.v1 as qs

job_tag = "test-4nt8"
traj = 1000

# Load gauge field and gauge transform
get_gf = qs.run_gf(job_tag, traj)
get_gt = qs.run_gt(job_tag, traj, get_gf)

# Generate selections
get_psel = qs.run_psel(job_tag, traj)
get_fsel = qs.run_fsel(job_tag, traj)
get_fselc = qs.run_fselc(job_tag, traj, get_fsel, get_psel)

# Generate wall-source info
get_wi = qs.run_wi(job_tag, traj)

# Compute eigenvectors
get_eig = qs.run_eig(job_tag, traj, get_gf)

# Access parameters
total_site = qs.get_param(job_tag, "total_site")
print(f"Lattice size: {total_site}")

q.end_with_mpi()
```
