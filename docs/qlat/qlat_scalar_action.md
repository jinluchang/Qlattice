# `qlat.scalar_action` — Scalar Field Action for HMC

Source: `qlat/qlat/scalar_action.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`ScalarAction` Class](#scalaraction-class)
   - [Constructor](#constructor)
   - [Parameter Accessors](#parameter-accessors)
   - [Action and Energy Methods](#action-and-energy-methods)
   - [HMC Molecular-Dynamics Methods](#hmc-molecular-dynamics-methods)
   - [Field Conversion Methods](#field-conversion-methods)
   - [Axial Current and Polar Field](#axial-current-and-polar-field)
3. [Examples](#examples)

---

## Overview

`scalar_action` implements a scalar field action on the lattice, parameterised
by a mass squared (`m_sq`) and two self-coupling constants (`lmbd`, `alpha`).
It is used by the HMC integrator in `qlat.hmc` to evolve scalar fields
through molecular dynamics.

The class wraps a C++ object (accessed through the `qlat.c` bindings) and
exposes all operations needed during an HMC trajectory: computing the action
on a single node, evaluating the Hamiltonian, setting forces, evolving fields,
and generating random momenta.

---

## `ScalarAction` Class

### Constructor

```python
ScalarAction(m_sq: float, lmbd: float, alpha: float)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `m_sq` | `float` | Mass squared (m^2) of the scalar field |
| `lmbd` | `float` | Quartic self-coupling (lambda) |
| `alpha` | `float` | Additional coupling parameter (alpha) |

```python
sa = q.ScalarAction(m_sq=0.1, lmbd=1.0, alpha=0.0)
```

---

### Parameter Accessors

#### `m_sq() -> float`

Return the mass-squared parameter.

#### `lmbd() -> float`

Return the quartic coupling parameter.

#### `alpha() -> float`

Return the additional coupling parameter (alpha).

---

### Action and Energy Methods

#### `action_node(sf: FieldBase) -> float`

Compute the scalar action contribution on the local MPI node for field `sf`.

#### `sum_sq(sf: FieldBase) -> float`

Return the sum of squared field values on the local node.

#### `hmc_m_hamilton_node(sf: FieldBase, masses: FieldBase) -> float`

Compute the Hamiltonian contribution from the kinetic (momentum) term on the
local node. `sf` is the momentum field and `masses` provides the
mode-dependent mass factors.

---

### HMC Molecular-Dynamics Methods

#### `hmc_set_force(sm_force: FieldBase, sf: FieldBase)`

Write the molecular-dynamics force derived from the action into `sm_force`,
evaluated at field configuration `sf`.

#### `hmc_field_evolve(sf_ft: FieldBase, sm_ft: FieldBase, masses: FieldBase, step_size: float)`

Advance the field (in Fourier space) by one leapfrog step:

- `sf_ft` — field in Fourier space (modified in place)
- `sm_ft` — conjugate momentum in Fourier space
- `masses` — mode-dependent mass factors
- `step_size` — integrator step size

#### `hmc_estimate_mass(masses: FieldBase, field_ft: FieldBase, force_ft: FieldBase, phi0: float)`

Estimate an effective mass for the Fourier-accelerated integrator based on the
current field and force.

#### `hmc_set_rand_momentum(sm_complex: FieldBase, masses: FieldBase, rs: RngState)`

Fill `sm_complex` with Gaussian random momenta drawn from `rs`, scaled by
`masses`.

#### `hmc_predict_field(field_ft: FieldBase, momentum_ft: FieldBase, masses: FieldBase, vev_sigma: float)`

Predict the evolved field from the current momentum and mass factors.

#### `to_mass_factor(sin_domega: FieldBase) -> FieldBase`

Convert the `sin_domega` field into a mass-factor field used by the
Fourier-accelerated integrator.

---

### Field Conversion Methods

#### `set_complex_from_double(cf: FieldBase, sf: FieldBase)`

Convert a real-valued field `sf` into a complex-valued field `cf`.

#### `set_double_from_complex(sf: FieldBase, cf: FieldBase)`

Convert a complex-valued field `cf` into a real-valued field `sf`.

---

### Axial Current and Polar Field

#### `axial_current_node(axial_current: FieldBase, sf: FieldBase)`

Compute the axial current on the local node from field `sf` and store the
result in `axial_current`.

#### `get_polar_field(polar_field: FieldBase, field: FieldBase)`

Decompose `field` into radial and angular (polar) components, storing the
result in `polar_field`.

---

## Examples

### Basic Construction and Parameter Query

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

sa = q.ScalarAction(m_sq=0.1, lmbd=1.0, alpha=0.0)
print(sa.m_sq())    # 0.1
print(sa.lmbd())    # 1.0

q.end_with_mpi()
```

### Computing the Action on a Field

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
sa = q.ScalarAction(m_sq=0.1, lmbd=1.0, alpha=0.0)
sf = q.Field(q.ElemTypeRealD, geo, 1)
sf.set_zero()

action = sa.action_node(sf)
print(f"action_node = {action}")

q.end_with_mpi()
```

### Copying an Action

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

sa1 = q.ScalarAction(0.1, 1.0, 0.0)
sa2 = q.ScalarAction(0.2, 2.0, 0.0)
sa2 @= sa1
print(sa2.m_sq())    # 0.1

q.end_with_mpi()
```
