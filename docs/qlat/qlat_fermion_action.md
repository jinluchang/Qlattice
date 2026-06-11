# `qlat.fermion_action` — Fermion Action Parameters for Domain-Wall Fermions

Source: `qlat/qlat/fermion_action.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`FermionAction` Class](#fermionaction-class)
   - [Constructor](#constructor)
   - [Attribute Accessors](#attribute-accessors)
   - [Assignment](#assignment)
3. [Examples](#examples)

---

## Overview

`fermion_action` provides the `FermionAction` class, which encapsulates the
parameters of a domain-wall fermion action. Two variants are supported:

- **Mobius** — a single `mobius_scale` parameter controls the extent of the
  fifth dimension. This is the standard formulation used in most lattice QCD
  simulations with domain-wall fermions.
- **ZMobius** — an explicit list of `omega` coefficients replaces the uniform
  Mobius scale, allowing finer control over the approximation to the
  sign function in the fifth dimension.

The class manages the lifetime of a C++ `FermionAction` object through the
`cdata` attribute and the `qlat.c` backend.

---

## `FermionAction` Class

### Constructor

```python
FermionAction(*, mass: float, ls: int, m5: float,
              mobius_scale: float = 1.0, omega: list = None)
```

Create a fermion action descriptor. If `omega` is `None`, a standard Mobius
action is created. If `omega` is provided (a list of length `ls`), a ZMobius
action is created instead.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `mass` | `float` | — | Quark mass |
| `ls` | `int` | — | Extent of the fifth dimension |
| `m5` | `float` | — | Domain-wall height (typically 1.0–1.8) |
| `mobius_scale` | `float` | `1.0` | Mobius scale parameter (Mobius variant only) |
| `omega` | `list[float]` or `None` | `None` | ZMobius coefficients; must have length `ls` |

**Raises:** `AssertionError` if parameter types are incorrect, or if
`len(omega) != ls` for the ZMobius variant.

---

### Attribute Accessors

#### `mass()`

```python
mass() -> float
```

Return the quark mass.

#### `ls()`

```python
ls() -> int
```

Return the extent of the fifth dimension.

#### `m5()`

```python
m5() -> float
```

Return the domain-wall height.

#### `omega()`

```python
omega() -> list[float]
```

Return the ZMobius omega coefficients. For a standard Mobius action this
returns the uniform omega values derived from `mobius_scale`.

#### `mobius_scale()`

```python
mobius_scale() -> float
```

Return the Mobius scale parameter.

---

### Assignment

#### `__imatmul__` (`@=`)

```python
fa1 @= fa2
```

Copy the contents of `fa2` into `fa1`. Both must be `FermionAction`
instances. This delegates to `qlat.c.set_fermion_action`.

---

## Examples

### Create a Mobius Fermion Action

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

fa = q.FermionAction(mass=0.01, ls=24, m5=1.8, mobius_scale=2.0)

print(f"mass = {fa.mass()}")
print(f"ls   = {fa.ls()}")
print(f"m5   = {fa.m5()}")
print(f"mobius_scale = {fa.mobius_scale()}")

q.end_with_mpi()
```

### Create a ZMobius Fermion Action

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

omega = [0.5 + 0.1 * i for i in range(12)]
fa = q.FermionAction(mass=0.01, ls=12, m5=1.8, omega=omega)

print(f"mass  = {fa.mass()}")
print(f"ls    = {fa.ls()}")
print(f"m5    = {fa.m5()}")
print(f"omega = {fa.omega()}")

q.end_with_mpi()
```

### Copy a Fermion Action

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

fa1 = q.FermionAction(mass=0.01, ls=24, m5=1.8, mobius_scale=2.0)
fa2 = q.FermionAction(mass=0.08, ls=16, m5=1.0, mobius_scale=1.0)

fa2 @= fa1
print(f"fa2 mass = {fa2.mass()}")  # 0.01
print(f"fa2 ls   = {fa2.ls()}")    # 24

q.end_with_mpi()
```
