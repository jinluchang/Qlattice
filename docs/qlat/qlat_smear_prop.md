# `qlat.smear_prop` — Spatial Propagator Smearing with MPI Communication

Source: `qlat/qlat/smear_prop.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Public API](#public-api)
   - [`prop_spatial_smear`](#prop_spatial_smear)
3. [Internal Helpers](#internal-helpers)
   - [`get_prop_spatial_smear_chunk_plan`](#get_prop_spatial_smear_chunk_plan)
   - [`prop_spatial_smear_chunk_planner`](#prop_spatial_smear_chunk_planner)
   - [`prop_spatial_smear_chunk`](#prop_spatial_smear_chunk)
4. [Examples](#examples)

---

## Overview

The `qlat.smear_prop` module implements **spatial Jacobi smearing** of
propagators and fermion fields in a distributed MPI environment. Unlike
`qlat.smear.prop_spatial_smear_no_comm`, which requires all time slices to
be local, this module handles the general case by redistributing fields
across nodes using `SelectedShufflePlan` so that each node can apply
`prop_spatial_smear_no_comm` to its local time slices.

Key design points:

- **Chunk-based processing** — Fields are processed in configurable chunks
  (default: 12 for `FermionField4d`, 1 for `Prop`) to control memory usage
  and plan reuse.
- **Plan caching** — The shuffle plan and redistributed gauge field are
  cached and reused when the gauge field and field count match, avoiding
  expensive re-planning on repeated calls.
- **Flop counting** — The main function reports floating-point operation
  counts for integration with the qlat timer infrastructure.
- **Transparent Prop support** — When a `Prop` (spin-color matrix field) is
  passed, it is automatically decomposed into 12 `FermionField4d` components,
  smeared, and reassembled.

The module depends on `qlat.c` for `GaugeField`, `Prop`, `FermionField4d`,
`PointsSelection`, `SelectedShufflePlan`, `prop_spatial_smear_no_comm`,
`mk_ff_list_from_prop`, and `mk_prop_from_ff_list`.

---

## Public API

### `prop_spatial_smear`

```python
prop_spatial_smear(
    ff_list: list[FermionField4d] | Prop | FermionField4d,
    gf: GaugeField,
    coef: float,
    step: int,
    mom: CoordinateD = None,
    *,
    chunk_size: int = None,
) -> list[FermionField4d] | Prop | FermionField4d
```

Perform spatial Jacobi smearing on `ff_list`. The return type matches the
input type:

| Input type | Return type |
|---|---|
| `list[FermionField4d]` | `list[FermionField4d]` |
| `list[Prop]` | `list[Prop]` |
| `Prop` | `Prop` |
| `FermionField4d` | `FermionField4d` |

The original `ff_list` is not modified; a new list (or object) is returned.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ff_list` | `list`, `Prop`, or `FermionField4d` | — | Fields to smear |
| `gf` | `GaugeField` | — | Gauge field (typically APE-smeared) |
| `coef` | `float` | — | Smearing coefficient (higher = less smearing) |
| `step` | `int` | — | Number of smearing iterations |
| `mom` | `CoordinateD` or `None` | `None` | Momentum for phase projection; `None` = zero |
| `chunk_size` | `int` or `None` | `None` | Chunk size for batching; default 12 for `FermionField4d`, 1 for `Prop` |
| **Returns** | `list/Prop/FermionField4d` | | Smeared fields (type matches input) |

When `step == 0`, a copy of the input is returned without smearing.

---

## Internal Helpers

### `get_prop_spatial_smear_chunk_plan`

```python
get_prop_spatial_smear_chunk_plan(gf: GaugeField, num_field: int) -> dict
```

Retrieve or compute a chunk smearing plan. Only one plan is cached at a time;
the cache is invalidated when the gauge field reference or field count changes.

| Parameter | Type | Description |
|---|---|---|
| `gf` | `GaugeField` | Gauge field |
| `num_field` | `int` | Number of fermion fields in the chunk |
| **Returns** | `dict` | Plan dictionary (see `prop_spatial_smear_chunk_planner`) |

---

### `prop_spatial_smear_chunk_planner`

```python
prop_spatial_smear_chunk_planner(gf: GaugeField, num_field: int) -> dict
```

Build a plan for chunk-based spatial smearing. The plan contains:

| Key | Type | Description |
|---|---|---|
| `gf` | `GaugeField` | Reference gauge field |
| `num_field` | `int` | Number of fields in the chunk |
| `s_gf_list` | `list` | Redistributed gauge field slices |
| `ssp` | `SelectedShufflePlan` | Shuffle plan for the fermion fields |
| `flops_per_step` | `int` | Flop count per smearing step |

---

### `prop_spatial_smear_chunk`

```python
prop_spatial_smear_chunk(
    ff_list: list[FermionField4d],
    plan: dict,
    coef: float,
    step: int,
    mom: CoordinateD,
) -> tuple[int, list[FermionField4d]]
```

Apply spatial smearing to a single chunk of fermion fields using a
precomputed plan. Fields are shuffled to per-time-slice layouts, smeared
in place with `prop_spatial_smear_no_comm`, and shuffled back.

| Parameter | Type | Description |
|---|---|---|
| `ff_list` | `list[FermionField4d]` | Chunk of fermion fields |
| `plan` | `dict` | Plan from `get_prop_spatial_smear_chunk_plan` |
| `coef` | `float` | Smearing coefficient |
| `step` | `int` | Number of smearing iterations |
| `mom` | `CoordinateD` | Momentum for phase projection |
| **Returns** | `tuple[int, list[FermionField4d]]` | `(flops, smeared_fields)` |

---

## Examples

### Basic Propagator Smearing

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Prepare APE-smeared gauge field
gf = q.GaugeField(geo)
gf.set_unit()

# Create a propagator and smear it
prop = q.Prop(geo)
prop.set_zero()

prop_smeared = q.smear_prop.prop_spatial_smear(
    prop, gf, coef=0.9375, step=54
)
print(f"Smeared propagator")

q.end_with_mpi()
```

### Smearing a List of FermionFields

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

gf = q.GaugeField(geo)
gf.set_unit()

# Create multiple fermion fields
ff_list = [q.FermionField4d(geo) for _ in range(6)]
for ff in ff_list:
    ff.set_zero()

# Smear with custom chunk size and momentum
mom = q.CoordinateD([0.1, 0.0, 0.0, 0.0])
smeared_list = q.smear_prop.prop_spatial_smear(
    ff_list, gf, coef=0.9375, step=30, mom=mom, chunk_size=6
)
print(f"Smeared {len(smeared_list)} fields")

q.end_with_mpi()
```

### Zero-Step Passthrough

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

gf = q.GaugeField(geo)
gf.set_unit()

ff = q.FermionField4d(geo)
ff.set_zero()

# step=0 returns a copy without smearing
ff_copy = q.smear_prop.prop_spatial_smear(ff, gf, coef=0.5, step=0)

q.end_with_mpi()
```
