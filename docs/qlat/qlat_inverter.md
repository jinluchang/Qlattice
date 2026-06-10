# `qlat.inverter` — Fermion-Matrix Inverter Framework

Source: `qlat/qlat/inverter.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`Inverter` Base Class](#inverter-base-class)
3. [`InverterDwfFreeField` Class](#inverterdwffreefield-class)
   - [Constructor](#constructor)
   - [Methods](#methods)
4. [`InverterDomainWall` Class](#inverterdomainwall-class)
   - [Constructor](#constructor-1)
   - [Methods](#methods-1)
   - [Solver Parameters](#solver-parameters)
5. [`InverterGaugeTransform` Class](#invertergaugetransform-class)
   - [Constructor](#constructor-2)
   - [Methods](#methods-2)
6. [`EigSystem` Class](#eigsystem-class)
7. [Module-Level Cache](#module-level-cache)
8. [Examples](#examples)

---

## Overview

`inverter` provides a lightweight framework for applying the inverse of the
domain-wall Dirac operator to `Prop` (propagator) fields. All concrete
inverters share a common interface: multiplying an inverter by a source
propagator returns the solution propagator.

The module defines:

- **`Inverter`** — empty base class establishing the interface convention.
- **`InverterDwfFreeField`** — analytic (free-field) inverse of the DWF
  operator; no gauge field is needed.
- **`InverterDomainWall`** — CG-based inversion backed by the C/C++ library;
  requires a gauge field and a `FermionAction`.
- **`InverterGaugeTransform`** — decorator that wraps any `Inverter`,
  applying a gauge transformation before and after inversion so that the
  result is in a different gauge.
- **`EigSystem`** — placeholder for future eigenvalue-system support.

All inverter classes support multiplication with a single `Prop` or a
`list` of `Prop` objects.

---

## `Inverter` Base Class

```python
class Inverter:
    pass
```

Empty base class. All concrete inverters inherit from it and implement
`__mul__` to apply the inverse Dirac operator.

---

## `InverterDwfFreeField` Class

Analytically inverts the free (no gauge field) domain-wall Dirac operator.
Useful for testing and as a reference solution.

### Constructor

```python
InverterDwfFreeField(*, mass, m5=1.0, momtwist=None, qtimer=TimerNone())
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `mass` | `float` | — | Fermion mass |
| `m5` | `float` | `1.0` | Fifth-dimensional mass parameter |
| `momtwist` | `CoordinateD` | `None` | Momentum twist (4-vector); defaults to zero |
| `qtimer` | `Timer` / `TimerNone` | `TimerNone()` | Optional timer for profiling |

### Methods

#### `__mul__(prop_src) -> Prop`

Apply the free-field inverse to a propagator source.

| Parameter | Type | Description |
|---|---|---|
| `prop_src` | `Prop` or `list[Prop]` | Source propagator(s) |

Returns the solution propagator, or a list of solution propagators if
`prop_src` is a list.

---

## `InverterDomainWall` Class

CG-based domain-wall fermion inverter backed by the C/C++ library. Requires
a gauge field and a `FermionAction` specification.

### Constructor

```python
InverterDomainWall(*, gf, fa, qtimer=TimerNone())
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Gauge field configuration |
| `fa` | `FermionAction` | — | Fermion action parameters (mass, Ls, m5, etc.) |
| `qtimer` | `Timer` / `TimerNone` | `TimerNone()` | Optional timer for profiling |

### Methods

#### `__mul__(prop_src) -> Prop`

Apply the domain-wall inverse to a propagator source using CG iteration.

| Parameter | Type | Description |
|---|---|---|
| `prop_src` | `Prop` or `list[Prop]` | Source propagator(s) |

Returns the solution propagator, or a list of solution propagators.

### Solver Parameters

#### `stop_rsd() -> float`

Return the current stopping residual for the CG solver.

#### `set_stop_rsd(stop_rsd) -> None`

Set the stopping residual for the CG solver.

#### `max_num_iter() -> int`

Return the maximum number of CG iterations.

#### `set_max_num_iter(max_num_iter) -> None`

Set the maximum number of CG iterations.

#### `max_mixed_precision_cycle() -> int`

Return the maximum number of mixed-precision refinement cycles.

#### `set_max_mixed_precision_cycle(max_mixed_precision_cycle) -> None`

Set the maximum number of mixed-precision refinement cycles.

---

## `InverterGaugeTransform` Class

Decorator inverter that wraps any `Inverter` with gauge transformations.
Given a gauge transform `gt`, it computes `gt_inv * src`, inverts using the
wrapped inverter, then applies `gt * sol` so the result is in the
transformed gauge.

### Constructor

```python
InverterGaugeTransform(*, inverter, gt, qtimer=TimerNone())
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inverter` | `Inverter` | — | Inner inverter to wrap |
| `gt` | `GaugeTransform` | — | Gauge transformation to apply |
| `qtimer` | `Timer` / `TimerNone` | `TimerNone()` | Optional timer for profiling |

### Methods

#### `__mul__(prop_src) -> Prop | FermionField4d`

Apply the gauge-transformed inversion: `gt * (inverter * (gt_inv * src))`.

| Parameter | Type | Description |
|---|---|---|
| `prop_src` | `Prop`, `FermionField4d`, or `list` | Source field(s) |

Returns the solution in the transformed gauge.

---

## `EigSystem` Class

```python
class EigSystem:
    pass
```

Placeholder class for future eigenvalue-system support. Currently empty.

---

## Module-Level Cache

```python
cache_inv = mk_cache("inv")
```

A module-level cache created via `mk_cache("inv")`. Used internally to
store and reuse intermediate inversion results.

---

## Examples

### Free-Field vs. Domain-Wall Inversion

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
rs = q.RngState("seed")

gf = q.GaugeField(geo)
gf.set_unit()

fa = q.FermionAction(mass=0.05, ls=16, m5=1.0)

qinv_free = q.InverterDwfFreeField(
    mass=fa.mass(), m5=fa.m5(), qtimer=q.Timer("InverterDwfFreeField")
)
qinv_dwf = q.InverterDomainWall(
    gf=gf, fa=fa, qtimer=q.Timer("InverterDomainWall")
)

src = q.Prop(geo)
src.set_rand(rs.split("src"))

sol_free = qinv_free * src
sol_dwf = qinv_dwf * src

sol_diff = sol_dwf.copy()
sol_diff -= sol_free
print(f"free  sol qnorm = {sol_free.qnorm():.6E}")
print(f"dwf   sol qnorm = {sol_dwf.qnorm():.6E}")
print(f"diff  sol qnorm = {sol_diff.qnorm():.3E}")

q.end_with_mpi()
```

### Gauge-Transformed Inversion

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
rs = q.RngState("gfix-test")

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf"), 0.05, 2)

fa = q.FermionAction(mass=0.1, ls=8, m5=1.0)
gt = q.GaugeTransform(geo)
gt.set_rand(rs.split("gt"), 0.1, 1)

gf_gfix = gt * gf

inv = q.InverterDomainWall(gf=gf_gfix, fa=fa, qtimer=q.Timer("inv"))
inv_gt = q.InverterGaugeTransform(
    inverter=q.InverterDomainWall(gf=gf, fa=fa, qtimer=q.Timer("inv-inner")),
    gt=gt,
    qtimer=q.Timer("inv-gt"),
)

src = q.Prop(geo)
src.set_rand(rs.split("src"))

sol = inv * src
sol_gt = inv_gt * src

print(f"sol    qnorm = {sol.qnorm():.6E}")
print(f"sol_gt qnorm = {sol_gt.qnorm():.6E}")

q.end_with_mpi()
```

### Inverting a List of Sources

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
rs = q.RngState("list-test")

gf = q.GaugeField(geo)
gf.set_unit()

fa = q.FermionAction(mass=0.05, ls=16, m5=1.0)
qinv = q.InverterDomainWall(gf=gf, fa=fa)

sources = []
for i in range(3):
    src = q.Prop(geo)
    src.set_rand(rs.split(f"src-{i}"))
    sources.append(src)

solutions = qinv * sources
for i, sol in enumerate(solutions):
    print(f"solution[{i}] qnorm = {sol.qnorm():.6E}")

q.end_with_mpi()
```
