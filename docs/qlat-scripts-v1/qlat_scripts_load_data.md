# `qlat_scripts.v1.load_data` — Propagator Data Loading and Cache Management

Source: `qlat/qlat_scripts/v1/load_data.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [AMA Value Construction](#ama-value-construction)
3. [Wall-Source Propagator Accessors](#wall-source-propagator-accessors)
4. [Point-Source Propagator Accessors](#point-source-propagator-accessors)
5. [Random U(1) Propagator Accessors](#random-u1-propagator-accessors)
6. [Flavor Mapping](#flavor-mapping)
7. [Generic Element Access](#generic-element-access)
8. [Cache Population](#cache-population)
9. [Propagator Lookup](#propagator-lookup)
10. [Data Loading Functions](#data-loading-functions)
11. [Gauge Link Access](#gauge-link-access)
12. [Top-Level Pipeline](#top-level-pipeline)
13. [Examples](#examples)

---

## Overview

This module provides the data-loading layer for auto-contractor evaluations in `qlat_scripts.v1`. It reads saved propagator data from disk (produced by `gen_data.py`), populates in-memory caches, and exposes a unified `get_prop(flavor, p_snk, p_src)` interface that returns AMA-corrected propagator values.

Key responsibilities:

- Load wall-source, point-source, and random U(1) propagators from field/QAR archives
- Construct AMA (all-mode-accumulation) values from multi-accuracy data
- Build lookup caches indexed by `(flavor, position, source_type, sink_type)`
- Handle gauge transformations and boundary condition corrections
- Provide a top-level `run_get_prop` that orchestrates all loading

---

## AMA Value Construction

### `get_prop_wsrc(prop_cache, inv_type, t_src, tag_snk_type)`

Constructs an AMA value for a wall-source propagator. Loads sloppy (`inv_acc=1`) and exact (`inv_acc=2`) data, combining them via `mk_ama_val` with the appropriate probability weights.

**Parameters:**
- `tag_snk_type` — one of `"wsrc_wsnk ; psel_ts"`, `"wsrc ; fsel"`, `"wsrc ; psel"`

### `get_prop_psrc(prop_cache, inv_type, xg_src, tag_snk_type)`

Constructs an AMA value for a point-source propagator. Combines up to four accuracy levels (0, 1, 2) via `mk_ama_val`. Handles the case where strange quark sloppy inversions may have fewer samples than light quark.

**Parameters:**
- `tag_snk_type` — one of `"psrc_wsnk ; psel_ts"`, `"psrc ; fsel"`, `"psrc ; psel"`

---

## Wall-Source Propagator Accessors

### `get_prop_wsnk_wsrc(prop_cache, inv_type, t_snk, t_src)`

Returns the wall-sink/wall-source propagator element at `(t_snk, t_src)`.

### `get_prop_psnk_wsrc_fsel(prop_cache, inv_type, xg_snk, t_src, fsel_pos_dict)`

Returns the point-sink/wall-source propagator element on the field selection.

### `get_prop_psnk_wsrc_psel(prop_cache, inv_type, xg_snk, t_src, psel_pos_dict)`

Returns the point-sink/wall-source propagator element on the point selection.

---

## Point-Source Propagator Accessors

### `get_prop_wsnk_psrc(prop_cache, inv_type, t_snk, xg_src)`

Returns the wall-sink/point-source propagator element.

### `get_prop_psnk_psrc_fsel(prop_cache, inv_type, xg_snk, xg_src, fsel_pos_dict)`

Returns the point-sink/point-source propagator element on the field selection.

### `get_prop_psnk_psrc_psel(prop_cache, inv_type, xg_snk, xg_src, psel_pos_dict)`

Returns the point-sink/point-source propagator element on the point selection.

---

## Random U(1) Propagator Accessors

### `get_prop_rand_u1_fsel(prop_cache, inv_type)`

Returns the averaged random U(1) propagator on the field selection (exact accuracy only).

### `get_prop_psnk_rand_u1_fsel(prop_cache, inv_type, xg_snk, xg_src, fsel_pos_dict)`

Returns a single element of the random U(1) propagator. Requires `xg_snk == xg_src`.

---

## Flavor Mapping

### `dict_flavor_inv_type`

Maps flavor string to inversion type integer:

| Flavor | `inv_type` |
|--------|-----------|
| `"l"`, `"u"`, `"d"` | 0 |
| `"s"` | 1 |
| `"c"` | 2 |

---

## Generic Element Access

### `f_get_elem_wm(field, pos_snk)`

Extracts a Wilson matrix element from a field, handling the zero-field (`int(0)`) case.

### `mk_get_elem_wm(field, pos_dict=None)`

Returns a `get(pos_snk)` function for extracting Wilson matrix elements. Supports both plain fields and `AmaVal` objects. If `pos_dict` is provided, translates coordinate keys to integer indices.

### `mk_field_norm_sqrt(field)`

Computes `sqrt(qnorm_field(field))` for norm-based element access.

### `mk_get_elem_norm(field, pos_dict=None)`

Returns a `get(pos_snk)` function for extracting propagator norms (square root). Used for importance-sampling weight calculations.

---

## Cache Population

### `populate_prop_idx_cache_wsrc_psel(...)`

Populates `prop_lookup_cache` and `prop_norm_lookup_cache` for wall-source propagators on the point selection.

### `populate_prop_idx_cache_wsrc_fsel(...)`

Populates caches for wall-source propagators on the field selection.

### `populate_prop_idx_cache_psrc_psel(...)`

Populates caches for point-source propagators on the point selection.

### `populate_prop_idx_cache_psrc_fsel(...)`

Populates caches for point-source propagators on the field selection.

### `populate_prop_idx_cache_rand_u1_fsel(...)`

Populates caches for random U(1) propagators on the field selection.

Each function creates cache entries keyed by `(flavor, pos_src, type_src, type_snk)`.

---

## Propagator Lookup

### `get_prop_lookup_snk_src(prop_lookup_cache, flavor, p_snk, p_src)`

Looks up a propagator element. Falls back to `g5_herm` conjugation if the reverse direction is cached.

**Parameters:**
- `p_snk`, `p_src` — tuples like `("point", xg_tuple)` or `("wall", t_int)`

### `get_prop_norm_lookup_snk_src(prop_norm_lookup_cache, flavor, p_snk, p_src)`

Same as above but returns the propagator norm (square root of qnorm).

---

## Data Loading Functions

### `load_prop_wsrc_psel(job_tag, traj, flavor, *, psel, fsel, gt)`

Loads wall-source point-selected propagators. Applies inverse gauge transformation. Handles `48I` strange quark boundary condition flip for MIRA data.

### `load_prop_wsrc_fsel(job_tag, traj, flavor, *, psel, fsel, gt)`

Loads wall-source field-selected propagators. Validates consistency with psel data.

### `load_prop_psrc_psel(job_tag, traj, flavor, *, psel, fsel)`

Loads point-source point-selected propagators. Computes AMA probabilities from actual sample counts.

### `load_prop_psrc_fsel(job_tag, traj, flavor, *, psel, fsel)`

Loads point-source field-selected propagators. Validates against psel data.

### `load_prop_rand_u1_fsel(job_tag, traj, flavor, *, psel, fsel)`

Loads random U(1) propagators and performs AMA reweighting across accuracy levels. Averages over all `n_rand_u1_fsel` sources.

### `load_gauge_hyp(job_tag, traj, *, gf_hyp)`

Loads HYP-smeared gauge field and its dagger into the gauge cache with expansion.

---

## Gauge Link Access

### `get_gauge_link_lookup_p_mu(prop_cache, tag, p, mu)`

Returns a `q.ColorMatrix` for the gauge link at position `p` in direction `mu`. Supports negative `mu` for backward links (adjoint).

---

## Top-Level Pipeline

### `run_get_prop(job_tag, traj, *, get_gf=None, get_gf_hyp=None, get_gt=None, get_psel, get_fsel, get_psel_smear=None, prop_types=None)`

Main entry point. Returns a lazy callable `get_prop(flavor, *args)` that:

1. Loads all required data on first invocation (gauge fields, propagators, gauge links)
2. Builds position dictionaries for psel and fsel
3. Populates lookup caches

**Default `prop_types`** (load order):
```
"wsrc psel s", "wsrc psel l", "wsrc fsel s", "wsrc fsel l",
"psrc psel s", "psrc psel l", "psrc fsel s", "psrc fsel l",
"rand_u1 fsel c", "rand_u1 fsel s", "rand_u1 fsel l", "gf hyp"
```

**`get_prop` interface:**
- `get_prop(flavor, p_snk, p_src)` — returns AMA propagator value
- `get_prop(flavor, p_snk, p_src, is_norm_sqrt=True)` — returns propagator norm
- `get_prop("U", tag, p, mu)` — returns gauge link matrix

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
get_psel = qs.run_psel(job_tag, traj)
get_fsel = qs.run_fsel(job_tag, traj)

get_prop = qs.run_get_prop(
    job_tag, traj,
    get_gf=get_gf, get_gt=get_gt,
    get_psel=get_psel, get_fsel=get_fsel,
)

# Lookup a wall-source to point-sink propagator
# p_src = ("wall", t_src)
# p_snk = ("point", xg_snk_tuple)
# result = get_prop("l", p_snk, p_src)

q.end_with_mpi()
```
