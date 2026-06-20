# `qlat.qlat_gpt` — GPT–Qlattice Interoperability

Source: `qlat/qlat_gpt.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Initialization](#initialization)
3. [Field Conversion](#field-conversion)
   - [Generic Converters](#generic-converters)
   - [Gauge Fields](#gauge-fields)
   - [Gauge Transforms](#gauge-transforms)
   - [Propagators](#propagators)
   - [Fermion Fields](#fermion-fields)
   - [Complex Fields](#complex-fields)
4. [Copy Plans](#copy-plans)
5. [Type Predicates](#type-predicates)
6. [Inverter and Eigensystem](#inverter-and-eigensystem)
   - [InverterGPT](#invertergpt)
   - [EigSystemGPT](#eigsystemgpt)
   - [EigSystemCompressedGPT](#eigsystemcompressedgpt)
7. [Gauge Fixing](#gauge-fixing)
   - [gauge_fix_coulomb](#gauge_fix_coulomb)
   - [check_gauge_fix_coulomb](#check_gauge_fix_coulomb)
   - [non_linear_cg](#non_linear_cg)
8. [Utilities](#utilities)
9. [Examples](#examples)

---

## Overview

`qlat_gpt` bridges the [GPT](https://github.com/lehner/gpt) (Grid) lattice
framework and Qlattice's field types. It provides:

- **MPI initialization** that coordinates GPT's Grid layout with Qlattice's
  node geometry.
- **Zero-copy-style field conversion** between GPT lattice objects and Qlattice
  `Field` / `GaugeField` / `Prop` / `FermionField4d` / `FieldComplexD` via
  cached `g.copy_plan` objects.
- **Inverter and eigensystem wrappers** (`InverterGPT`, `EigSystemGPT`,
  `EigSystemCompressedGPT`) that adapt GPT's linear-algebra solvers to
  Qlattice's `Inverter` / `EigSystem` interfaces.
- **Coulomb gauge fixing** (`gauge_fix_coulomb`) built on GPT's Landau fixing
  with split-grid acceleration and a custom non-linear CG optimizer.
- **NERSC gauge field I/O** via GPT's format support.

---

## Initialization

### `begin_with_gpt()`

Initialize both GPT and Qlattice runtimes. Must be called **before** any other
Qlattice or GPT operation.

1. Creates a GPT grid (from `--grid` argument or a default 288×288×288×576
   lattice).
2. Derives `size_node`, `coor_node`, and `id_node` from the GPT grid's MPI
   layout.
3. Calls `q.begin()` and `q.set_comm()` with the appropriate MPI communicator.

### `end_with_gpt()`

Tear down both runtimes. Frees the Qlattice communicator and calls `q.end()`.

```python
import qlat_gpt as qg
qg.begin_with_gpt()
# ... work ...
qg.end_with_gpt()
```

---

## Field Conversion

### Generic Converters

#### `qlat_from_gpt(gpt_obj) → object`

Auto-detect the GPT object type and convert it to the corresponding Qlattice
type:

| GPT type | Qlattice type |
|---|---|
| `list[4]` of `mcolor` | `GaugeField` |
| `mcolor` (single) | `GaugeTransform` |
| `mspincolor` | `Prop` |
| `vspincolor` | `FermionField4d` |
| `list` of `complex` | `FieldComplexD` |
| `list` of any above | `list` of converted elements |

Raises `Exception` for unrecognized types.

#### `gpt_from_qlat(obj) → object`

Reverse of `qlat_from_gpt`. Auto-detects Qlattice types (`Prop`,
`GaugeTransform`, `GaugeField`, `FermionField4d`, `FieldComplexD`) and converts
to the GPT equivalent. Lists are converted element-wise.

### Gauge Fields

#### `qlat_from_gpt_gauge_field(gpt_gf) → GaugeField`

Convert a list of 4 GPT `mcolor` fields (one per direction μ=0,1,2,3) into a
single Qlattice `GaugeField` with multiplicity 4. Each direction is converted
individually via a cached copy plan, then merged with `q.merge_fields`.

#### `gpt_from_qlat_gauge_field(gf) → list[mcolor]`

Split a Qlattice `GaugeField` into 4 single-direction `FieldColorMatrix`
objects, then convert each to a GPT `mcolor` field.

### Gauge Transforms

#### `qlat_from_gpt_gauge_transform(gpt_gt) → GaugeTransform`

Convert a single GPT `mcolor` field into a Qlattice `GaugeTransform`.

#### `gpt_from_qlat_gauge_transform(gt) → mcolor`

Convert a Qlattice `GaugeTransform` into a GPT `mcolor` field.

### Propagators

#### `qlat_from_gpt_prop(gpt_prop) → Prop`

Convert a GPT `mspincolor` propagator to a Qlattice `Prop` (WilsonMatrix).
Internally converts through `mspincolor → WilsonMatrix` via
`q.convert_wm_from_mspincolor`.

#### `gpt_from_qlat_prop(prop_wm) → mspincolor`

Convert a Qlattice `Prop` to a GPT `mspincolor`. Uses
`q.convert_mspincolor_from_wm` before the copy plan.

### Fermion Fields

#### `qlat_from_gpt_ff4d(gpt_ff) → FermionField4d`

Convert a GPT `vspincolor` field to a Qlattice `FermionField4d`.

#### `gpt_from_qlat_ff4d(ff) → vspincolor`

Convert a Qlattice `FermionField4d` to a GPT `vspincolor`.

### Complex Fields

#### `qlat_from_gpt_complex(gpt_fcs) → FieldComplexD`

Convert a list of GPT `complex` fields into a single Qlattice `FieldComplexD`.
If the list has one element, returns a multiplicity-1 field. Otherwise, converts
each element separately and merges into a higher-multiplicity field.

#### `gpt_from_qlat_complex(fc) → list[complex]`

Split a Qlattice `FieldComplexD` (possibly multi-multiplicity) into individual
GPT `complex` fields.

---

## Copy Plans

The conversion functions all use GPT's `g.copy_plan` mechanism for efficient
data transfer. Plans are created once and cached.

#### `mk_qlat_gpt_copy_plan_key(ctype, total_site, multiplicity, tag) → str`

Build a cache key string from element type name, total site, multiplicity, and
direction tag (`"qlat_from_gpt"` or `"gpt_from_qlat"`).

#### `mk_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag) → copy_plan`

Create a GPT `copy_plan` that maps between a GPT lattice field (in
lexicographic order) and a Qlattice `memoryview` buffer. The plan is built by:

1. Creating temporary GPT and Qlattice fields.
2. Collecting lexicographic coordinates from GPT.
3. Setting up `global_memory_view` on the Qlattice buffer side.
4. Associating source/destination views.
5. Finalizing with `local_only=True`.

#### `get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag) → copy_plan`

Retrieve a cached plan or create and cache a new one.

#### `mk_gpt_field(ctype, geo) → lattice`

Create a GPT lattice field matching the given Qlattice element type:

| `ctype` | GPT type |
|---|---|
| `ElemTypeColorMatrix` | `g.mcolor(grid)` |
| `ElemTypeWilsonMatrix` | `g.mspincolor(grid)` |
| `ElemTypeWilsonVector` | `g.vspincolor(grid)` |
| `ElemTypeComplexD` | `g.complex(grid)` |

---

## Type Predicates

Functions to identify GPT object types by inspecting `describe()` strings.

| Function | Matches |
|---|---|
| `is_gpt_prop(obj)` | `ot_matrix_spin_color(4,3);none` |
| `is_gpt_ff4d(obj)` | `ot_vector_spin_color(4,3);none` |
| `is_gpt_gauge_field(obj)` | `list[4]` of `ot_matrix_su_n_fundamental_group(3);none` |
| `is_gpt_gauge_transform(obj)` | `ot_matrix_su_n_fundamental_group(3);none` (single) |
| `is_gpt_complex(obj)` | `list` of `ot_complex_additive_group;none` |

---

## Inverter and Eigensystem

### InverterGPT

```python
class InverterGPT(q.Inverter):
    def __init__(self, *, inverter, qtimer=..., gpt_qtimer=...)
```

Wraps a GPT inverter (e.g. `g.algorithms.inverter`) as a Qlattice `Inverter`.
Supports the `inv * prop_src` syntax:

1. Converts the Qlattice source to GPT format via `gpt_from_qlat`.
2. Applies the GPT inverter: `g.eval(inverter * src)`.
3. Converts the result back via `qlat_from_gpt`.

Accepts `Prop`, `FermionField4d`, or `list` as the source operand.

**Parameters:**

| Name | Type | Description |
|---|---|---|
| `inverter` | GPT inverter | The GPT linear solver object |
| `qtimer` | `q.Timer` / `q.TimerNone` | Timer for the full conversion+solve cycle |
| `gpt_qtimer` | `q.Timer` / `q.TimerNone` | Timer for the GPT solve step only |

### EigSystemGPT

```python
class EigSystemGPT(q.EigSystem):
    def __init__(self, *, evec=None, evals=None)
```

Stores and I/Os a full eigensystem (eigenvectors + eigenvalues) using GPT's
serialization.

**Methods:**

| Method | Description |
|---|---|
| `load(path)` | Load eigenvectors and eigenvalues from a GPT file |
| `save(path)` | Save to a GPT file (no-op if `path` is `None`) |

### EigSystemCompressedGPT

```python
class EigSystemCompressedGPT(q.EigSystem):
    def __init__(self, *, basis=None, cevec=None, evals=None)
```

Stores a compressed eigensystem (basis + compressed eigenvectors + eigenvalues).

**Methods:**

| Method | Description |
|---|---|
| `load(path, *, total_site, fermion_params)` | Load with grid reconstruction from fermion parameters |
| `save(path, *, nsingle, mpi)` | Save with `g.format.cevec` compression |

`load` requires `total_site` and `fermion_params` to reconstruct the fermion
grid (calls `get_fgrid`). `save` takes `nsingle` (number of single-precision
components) and `mpi` (MPI layout for the format).

---

## Gauge Fixing

### `gauge_fix_coulomb`

```python
def gauge_fix_coulomb(
    gf,
    *,
    gt=None,
    mpi_split=None,
    maxiter_gd=10,
    maxiter_cg=200,
    maxcycle_cg=50000,
    log_every=1,
    eps=1e-12,
    step=0.3,
    step_gd=0.1,
    rng_seed=None,
) → GaugeTransform
```

Fix a gauge field to Coulomb gauge using split-grid optimization.

**Algorithm:**

1. Convert the Qlattice `GaugeField` to GPT and separate into time slices.
2. Split time slices across MPI ranks using `g.split` with configurable
   `mpi_split`.
3. For each local time slice, run a two-stage optimizer:
   - **Gradient descent** (`maxiter_gd` iterations) for initial convergence.
   - **Non-linear CG** (`maxiter_cg` iterations per cycle, up to
     `maxcycle_cg` cycles) for refinement.
4. Use Fourier-accelerated gradients (`inverse_phat_square`).
5. Unsplit, project to SU(3), and verify convergence.
6. Return the gauge transformation as a Qlattice `GaugeTransform`.

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field |
| `gt` | `GaugeTransform` | `None` | Initial guess (unit if `None`) |
| `mpi_split` | `list[int]` | `[1,1,1]` | Spatial MPI split for the sub-grid |
| `maxiter_gd` | `int` | `10` | Max gradient descent iterations |
| `maxiter_cg` | `int` | `200` | Max CG iterations per cycle |
| `maxcycle_cg` | `int` | `50000` | Max CG restart cycles |
| `log_every` | `int` | `1` | Log functional every N iterations |
| `eps` | `float` | `1e-12` | Convergence threshold |
| `step` | `float` | `0.3` | CG step size |
| `step_gd` | `float` | `0.1` | Gradient descent step size |
| `rng_seed` | `str` | `None` | Random seed (no randomization if `None`) |

### `check_gauge_fix_coulomb`

```python
def check_gauge_fix_coulomb(gf, gt, eps=1e-12) → bool
```

Verify that a gauge transformation `gt` fixes `gf` to Coulomb gauge. Computes
the average θ (gradient norm squared per site per color) over all time slices.
Returns `True` if `θ < eps`.

### `non_linear_cg`

```python
class non_linear_cg(g.algorithms.base_iterative):
    def __init__(self, params)
```

Non-linear conjugate gradient optimizer built on GPT's optimization framework.
Used internally by `gauge_fix_coulomb`.

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `eps` | `float` | `1e-8` | Convergence threshold on |df|/sqrt(dof) |
| `maxiter` | `int` | `1000` | Maximum iterations |
| `step` | `float` | `1e-3` | Base step size for line search |
| `log_functional_every` | `int` | `10` | Log interval |
| `line_search` | callable | `line_search_quadratic` | Line search function |
| `beta` | callable | `fletcher_reeves` | CG β formula (default Fletcher–Reeves; `gauge_fix_coulomb` uses Polak–Ribière) |
| `max_c` | `int` | `3` | Maximum line search step multiplier |

Returns a callable `opt(x, dx)` that minimizes `f` over group-valued `x`.

### `line_search_quadratic`

```python
def line_search_quadratic(s, x, dx, dv0, df, step, *, max_c=3)
```

Quadratic line search used by `non_linear_cg`. Fits `f(x) = a + b*(x-c)^2`
from directional derivative sign changes to estimate the minimizer along
direction `s`. Returns the step multiplier `c` (or `None` on failure).

---

## Utilities

### `mk_grid(geo=None) → g.grid`

Create a GPT grid object. If `geo` is provided, uses its `total_site`.
Otherwise reads `--grid` from command-line arguments, falling back to a default
288×288×288×576 lattice.

### `get_fgrid(total_site, fermion_params) → grid`

Reconstruct a GPT fermion grid (possibly even-odd) from Qlattice geometry and
GPT fermion parameters. Supports both `mobius` and `zmobius` (detected by the
presence of `"omega"` in `fermion_params`).

### `save_gauge_field(gf, path)`

Save a Qlattice `GaugeField` to disk in NERSC format via GPT.

### `load_gauge_field(path) → GaugeField`

Load a gauge field from a NERSC-format file via GPT and return as a Qlattice
`GaugeField`.

---

## Examples

### Basic setup and field conversion

```python
import qlat_gpt as qg
qg.begin_with_gpt()

import qlat as q
import gpt as g

# Load a gauge field from disk
gf = qg.load_gauge_field("config.nersc")

# Convert to GPT for analysis
gpt_gf = qg.gpt_from_qlat(gf)

# ... do GPT-level operations ...

# Convert back to Qlattice
gf2 = qg.qlat_from_gpt(gpt_gf)

qg.end_with_gpt()
```

### Using InverterGPT

```python
import qlat_gpt as qg
qg.begin_with_gpt()

import qlat as q
import gpt as g

# Set up a GPT inverter
gpt_gf = qg.gpt_from_qlat(gf)
q = g.qcd.fermion.mobius(gpt_gf, params)
inv = g.algorithms.inverter.cg({"eps": 1e-8, "maxiter": 1000})
prop_sol = g.eval(inv * q.propagator(prop_src))

# Or wrap as Qlattice inverter
qlat_inv = qg.InverterGPT(inverter=inv * q.propagator)
prop_src_qlat = q.Prop(gf.geo)
# ... set up source ...
prop_sol_qlat = qlat_inv * prop_src_qlat

qg.end_with_gpt()
```

### Coulomb gauge fixing

```python
import qlat_gpt as qg
qg.begin_with_gpt()

import qlat as q

gf = qg.load_gauge_field("config.nersc")

# Fix to Coulomb gauge
gt = qg.gauge_fix_coulomb(gf, maxiter_cg=200, eps=1e-12)

# Verify
assert qg.check_gauge_fix_coulomb(gf, gt, eps=1e-12)

qg.end_with_gpt()
```
