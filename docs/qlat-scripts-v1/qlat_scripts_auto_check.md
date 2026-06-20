# `qlat_scripts.v1.auto_check` — Checker Routines for Propagator Validation

Source: `qlat/qlat_scripts/v1/auto_check.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Wall-Source Checker](#wall-source-checker)
3. [Point-Source Checker](#point-source-checker)
4. [Propagator Loading](#propagator-loading)
5. [Propagator Lookup](#propagator-lookup)
6. [Examples](#examples)

---

## Overview

This module provides simplified checker routines for validating wall-source and point-source propagator computations. Unlike the full `gen_data` pipeline, these functions perform exact-accuracy inversions at all time slices or all lattice points without AMA multi-accuracy or importance sampling. They are intended for correctness checks and testing.

Key features:

- Wall-source propagator generation at every time slice with exact accuracy
- Point-source propagator generation at every lattice point with exact accuracy
- Propagator loading from disk into a cache
- A `get_prop(flavor, p_snk, p_src)` lookup function that returns Wilson matrix elements for all source/sink combinations (point-point, point-wall, wall-point, wall-wall)

---

## Wall-Source Checker

### `get_all_points(total_site)`

Returns a plain Python `list` of all lattice points as `q.Coordinate` objects.

### `get_all_points_psel(total_site)`

Returns a `q.PointsSelection` containing every point on the lattice (full selection).

### `run_get_inverter_checker(job_tag, traj, *, inv_type, get_gf, get_gt=None, get_eig=None)`

Pre-computes and caches the exact-accuracy (`inv_acc=2`) inverter. Used internally by the checker routines.

### `compute_prop_wsrc_checker(job_tag, tslice, inv_type, inv_acc, *, idx, gf, gt, sfw, path_sp, eig, finished_tags)`

Performs a single wall-source inversion at `tslice`. Skips if `tag` is already in `finished_tags`. Saves the selected-field propagator and wall-sink projection.

### `compute_prop_wsrc_all_checker(job_tag, traj, *, inv_type, gf, gt, eig)`

Runs wall-source inversions for all time slices. Saves outputs to:
- `{job_tag}/prop-wsrc-{flavor}/traj-{traj}/` (selected-field data)
- `{job_tag}/psel-prop-wsrc-{flavor}/traj-{traj}/` (point-selected data)

### `run_prop_wsrc_checker(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt)`

Top-level entry point for wall-source checker. Skips if output already exists. Uses file-based locking (`q.obtain_lock`) for mutual exclusion.

---

## Point-Source Checker

### `compute_prop_psrc_checker(job_tag, xg_src, inv_type, inv_acc, *, idx, gf, gt, sfw, eig, finished_tags)`

Performs a single point-source inversion at position `xg_src`. Skips if `tag` is already in `finished_tags`.

### `compute_prop_psrc_all_checker(job_tag, traj, *, inv_type, gf, gt, eig)`

Runs point-source inversions for all lattice points (from `get_all_points`). Saves output to `{job_tag}/prop-psrc-{flavor}/traj-{traj}/`.

### `run_prop_psrc_checker(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt)`

Top-level entry point for point-source checker. Skips if output already exists.

---

## Propagator Loading

### `load_prop_psrc(job_tag, traj, inv_type)`

Loads all point-source propagators from disk into the `prop_cache` under key `"psnk-psrc"`. Each propagator is projected onto the full point selection via `q.PselProp`.

### `load_prop_wsrc(job_tag, traj, inv_type)`

Loads all wall-source propagators from disk into the `prop_cache` under keys `"psnk-wsrc"` (point-sink) and `"wsnk-wsrc"` (wall-sink).

### `run_get_prop_checker(job_tag, traj, *, get_gf, get_gt)`

Returns a lazy callable `get_prop(flavor, p_snk, p_src)` that performs all loading on first invocation. The lookup function supports:

| Source | Sink | Lookup |
|--------|------|--------|
| point  | point | `cache["psnk-psrc"][src_idx].get_elem_wm(snk_idx)` |
| point  | wall  | `cache["psnk-wsrc"][t_src].get_elem_wm(snk_idx)` |
| wall   | point | `g5_herm(cache["psnk-wsrc"][t_src].get_elem_wm(src_idx))` |
| wall   | wall  | `cache["wsnk-wsrc"][t_src].get_elem_wm(t_snk)` |

---

## Propagator Lookup

The `get_prop` function returned by `run_get_prop_checker` uses the convention:

- `p_snk` and `p_src` are `(type_tag, position)` tuples
- `type_tag` is `"point"` or `"wall"`
- For `"point"`, position is a `q.Coordinate`; for `"wall"`, position is a time-slice integer
- Flavor is `"l"` (light) or `"s"` (strange)

---

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])
import qlat_scripts.v1 as qs

job_tag = "test-4nt8"
traj = 1000

get_gf = qs.run_gf(job_tag, traj)
get_gt = qs.run_gt(job_tag, traj, get_gf)

# Run wall-source checker for light quarks
qs.run_prop_wsrc_checker(job_tag, traj, inv_type=0,
                         get_gf=get_gf, get_eig=None, get_gt=get_gt)

# Run point-source checker for strange quarks
qs.run_prop_psrc_checker(job_tag, traj, inv_type=1,
                         get_gf=get_gf, get_eig=None, get_gt=get_gt)

# Load and use propagators
get_prop = qs.run_get_prop_checker(job_tag, traj, get_gf=get_gf, get_gt=get_gt)
# get_prop("l", ("point", xg_snk), ("wall", t_src)) after first call

q.end_with_mpi()
```
