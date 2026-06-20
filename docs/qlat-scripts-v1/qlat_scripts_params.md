# `qlat_scripts.v1.params` — Ensemble Simulation Parameters

Source: `qlat/qlat_scripts/v1/params.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Trajectory Lists](#trajectory-lists)
3. [Field Selection Parameters](#field-selection-parameters)
4. [Wall-Source AMA Parameters](#wall-source-ama-parameters)
5. [Smeared Selection Parameters](#smeared-selection-parameters)
6. [Point-Source AMA Parameters](#point-source-ama-parameters)
7. [Random U(1) Parameters](#random-u1-parameters)
8. [Propagator Smearing Parameters](#propagator-smearing-parameters)
9. [Gauge Smearing Parameters](#gauge-smearing-parameters)
10. [Fermion Parameters](#fermion-parameters)
11. [Meson Observables](#meson-observables)
12. [Examples](#examples)

---

## Overview

This module defines simulation parameters for each ensemble in `dict_params`. It extends the base parameters defined in `rbc_ukqcd_params` with values specific to the data production pipeline: trajectory lists, selection rates, AMA probabilities, smearing coefficients, and meson observable configurations.

Parameters are set via `dict_params[ensemble][tag]` (direct assignment) or `set_param(ensemble, tag, value)` (function call). All parameters are accessed at runtime via `get_param(job_tag, tag)`.

Supported ensembles include: `test-4nt8`, `test-4nt16`, `test-8nt16`, `24D`, `32D`, `32Dfine`, `24DH`, `48I`, `64I`, `64I-pq`, `16IH2`, `32IfineH`, `32IcoarseH1`, `24IH1`, `24IH2`, `24IH3`, `32IH1`, `32IH2`, `32IH3`.

---

## Trajectory Lists

**Tag:** `traj_list`

Each ensemble specifies a list of trajectory numbers available for analysis. These are Python `range` objects converted to lists.

| Ensemble | Trajectories |
|----------|-------------|
| `test-4nt8` | 1000, 1040, ..., 1360 |
| `48I` | 2300, 2299, ..., 501 |
| `64I` | 500, 510, ..., 5990 |
| `24D` | 5100, 5090, ..., 1010 |
| `32D` | 500, 510, ..., 2990 |

---

## Field Selection Parameters

### `field_selection_fsel_rate`

Fraction of spatial sites selected per time slice for field selection.

| Ensemble | Rate |
|----------|------|
| `test-*` | 1/16 |
| `24D` | 1/16 |
| `48I`, `64I` | 1/32 |

### `field_selection_psel_rate`

Fraction of total lattice volume for point selection. Used by `get_n_points_psel`.

### `field_selection_fsel_psrc_prop_norm_threshold`

Threshold for adaptive point-source field selection in `gen_data.py`. Sites with propagator norm above this value are probabilistically added to the field selection.

| Ensemble | Threshold |
|----------|-----------|
| `test-*` | 1e-3 |
| `24D` | 1e-4 |
| `48I` | 4e-5 |
| `64I` | 2e-5 |

### `n_points_psel`

Number of points in the point selection (alternative to `field_selection_psel_rate`).

---

## Wall-Source AMA Parameters

### `prob_exact_wsrc`

Probability that a given time slice receives an exact-accuracy wall-source inversion. For ensembles using `n_exact_wsrc` instead, the probability is computed as `1 - (1 - 1/T)^n`.

| Ensemble | Probability |
|----------|------------|
| `test-4nt16` | 1/8 |
| `48I` | 1/48 |
| `64I` | 1/64 |
| `24D` | 1/32 |

---

## Smeared Selection Parameters

### `n_per_tslice_smear`

Number of points per time slice in the smeared point selection.

### `prob_acc_1_smear` / `prob_acc_2_smear`

AMA probabilities for smeared-source medium and exact accuracy inversions.

---

## Point-Source AMA Parameters

### `prob_acc_1_psrc` / `prob_acc_2_psrc`

AMA probabilities for point-source medium and exact accuracy inversions.

| Ensemble | `prob_acc_1` | `prob_acc_2` |
|----------|-------------|-------------|
| `test-*` | 1/4 | 1/16 |
| `48I` | 1/32 | 1/128 |
| `64I` | 1/32 | 1/128 |

---

## Random U(1) Parameters

### `n_rand_u1_fsel`

Number of random U(1) volume sources per trajectory.

| Ensemble | Count |
|----------|-------|
| `test-*` | 4 |
| `16IH2` | 16 |
| Most others | 64 |

### `prob_acc_1_rand_u1` / `prob_acc_2_rand_u1`

AMA probabilities for random U(1) medium and exact accuracy inversions.

---

## Propagator Smearing Parameters

### `prop_smear_coef`

APE smearing coefficient for propagator source smearing.

| Ensemble | Coefficient |
|----------|------------|
| Most | 0.9375 |

### `prop_smear_step`

Number of APE smearing steps for propagator sources.

| Ensemble | Steps |
|----------|-------|
| `24D` | 10 |
| `48I` | 29 |
| `64I` | 54 |
| `32IfineH` | 96 |

---

## Gauge Smearing Parameters

### `gf_ape_smear_coef` / `gf_ape_smear_step`

APE smearing parameters for gauge field smearing (used by `run_gf_ape`). Coefficient 0.5, 30 steps for all ensembles.

---

## Fermion Parameters

### `fermion_params`

Nested dictionary of fermion action parameters indexed by `[inv_type][inv_acc]`. The `Ls` (domain-wall extent) parameter can be overridden per ensemble.

---

## Meson Observables

### `meson_tensor_tsep`

Source-sink separation for meson tensor measurements.

| Ensemble | Separation |
|----------|-----------|
| `test-4nt8` | 1 |
| `24D` | 8 |
| `48I` | 12 |
| `64I` | 18 |

### `meson_jwjj_threshold`

Threshold for JWJJ (current-current) operator selection.

### `meson_tsep_list`

List of source-sink separations for meson correlator measurements.

---

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])
import qlat_scripts.v1 as qs

# Access parameters
job_tag = "test-4nt8"
total_site = qs.get_param(job_tag, "total_site")
print(f"Lattice size: {total_site}")

traj_list = qs.get_param(job_tag, "traj_list")
print(f"Trajectories: {traj_list[:5]}...")

fsel_rate = qs.get_param(job_tag, "field_selection_fsel_rate")
print(f"Field selection rate: {fsel_rate}")

q.end_with_mpi()
```
