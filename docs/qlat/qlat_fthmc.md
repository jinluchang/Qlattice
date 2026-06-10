# `qlat.fthmc` — Flow-Time Hamiltonian Monte Carlo Utilities

Source: `qlat/qlat/fthmc.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`FlowInfo` Class](#flowinfo-class)
   - [Constructor](#constructor)
   - [Methods](#methods)
3. [Free Functions](#free-functions)
   - [Flow Evolution](#flow-evolution)
   - [Flowed Hamiltonian](#flowed-hamiltonian)
   - [Flowed Force](#flowed-force)
4. [Examples](#examples)

---

## Overview

`fthmc` provides building blocks for flow-time Hamiltonian Monte Carlo
(ftHMC) simulations. In ftHMC the gauge action used in the Metropolis step
is evaluated on a *gradient-flowed* gauge field rather than the original
one, which can suppress short-distance fluctuations and improve
acceptance rates.

Key components:

- **`FlowInfo`** — accumulates a sequence of gradient-flow steps
  (direction, step size, number of steps) to be applied during HMC.
- **`gf_flow` / `gf_flow_inv`** — apply or invert the flow defined by a
  `FlowInfo` object.
- **`gf_hamilton_flowed_node`** — evaluate the gauge Hamiltonian on the
  flowed field (node-local).
- **`set_gm_force_flowed` / `set_gm_force_flowed_no_det`** — compute the
  molecular-dynamics force including the Jacobian of the flow.

All functions delegate to the C++ backend via `qlat.c`.

---

## `FlowInfo` Class

`FlowInfo` describes a sequence of gradient-flow operations to be applied
to a gauge field. Each call to `add_flow` or `add_rand_order_flow` appends
a step to the sequence.

### Constructor

```python
FlowInfo()
```

Create an empty flow descriptor with no steps.

### Methods

#### `add_flow(eo, mu, epsilon, flow_size=1)`

```python
add_flow(eo: int, mu: int, epsilon: float, flow_size: int = 1) -> None
```

Append a deterministic flow step.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `eo` | `int` | — | Even/odd parity selector |
| `mu` | `int` | — | Lattice direction (0–3) |
| `epsilon` | `float` | — | Flow step size |
| `flow_size` | `int` | `1` | Number of sub-steps |

#### `add_rand_order_flow(rng, epsilon, *args)`

```python
add_rand_order_flow(rng: RngState, epsilon: float, epsilon2: float = None) -> None
```

Append a randomised-order flow step. If `epsilon2` is provided, a
two-epsilon variant is used.

| Parameter | Type | Description |
|---|---|---|
| `rng` | `RngState` | Random state for ordering |
| `epsilon` | `float` | Flow step size |
| `epsilon2` | `float` (optional) | Second step size for the two-epsilon variant |

#### `show()`

```python
show() -> str
```

Return a human-readable string describing the accumulated flow steps.

---

## Free Functions

### Flow Evolution

#### `gf_flow`

```python
gf_flow(gf: GaugeField, gf0: GaugeField, fi: FlowInfo) -> None
```

Apply the flow described by `fi` to `gf0` and store the result in `gf`.
Both `gf` and `gf0` must be `GaugeField` instances.

#### `gf_flow_inv`

```python
gf_flow_inv(gf: GaugeField, gf1: GaugeField, fi: FlowInfo) -> None
```

Apply the *inverse* of the flow described by `fi` to `gf1` and store the
result in `gf`. This reverses the gradient-flow evolution.

---

### Flowed Hamiltonian

#### `gf_hamilton_flowed_node`

```python
gf_hamilton_flowed_node(gf0: GaugeField, ga: GaugeAction, fi: FlowInfo) -> float
```

Evaluate the gauge Hamiltonian (using action `ga`) on the flowed
counterpart of `gf0`. The flow is defined by `fi`. Returns the
node-local (non-communicating) Hamiltonian contribution.

---

### Flowed Force

#### `set_gm_force_flowed`

```python
set_gm_force_flowed(
    gm_force: GaugeMomentum,
    gf0: GaugeField,
    ga: GaugeAction,
    fi: FlowInfo,
) -> None
```

Compute the molecular-dynamics force on the gauge field including the
Jacobian of the gradient flow, and store it in `gm_force`.

#### `set_gm_force_flowed_no_det`

```python
set_gm_force_flowed_no_det(
    gm_force: GaugeMomentum,
    gm_force_pre: GaugeMomentum,
    gf0: GaugeField,
    fi: FlowInfo,
) -> None
```

Compute the flowed force *without* the determinant Jacobian
contribution. `gm_force_pre` provides a pre-computed force that is
transformed by the flow Jacobian and added to `gm_force`.

---

## Examples

### Basic Flow and Inverse

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf0 = q.GaugeField(geo)
gf0.unitarize()

fi = q.FlowInfo()
fi.add_flow(0, 0, 0.1, flow_size=1)
fi.add_flow(0, 1, 0.1, flow_size=1)

gf1 = q.GaugeField(geo)
q.gf_flow(gf1, gf0, fi)

gf2 = q.GaugeField(geo)
q.gf_flow_inv(gf2, gf1, fi)

# gf2 should be close to gf0
print("Flow and inverse flow completed.")

q.end_with_mpi()
```

### Flowed Hamiltonian

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf0 = q.GaugeField(geo)
gf0.unitarize()

ga = q.GaugeAction(6.0)

fi = q.FlowInfo()
fi.add_flow(0, 0, 0.05, flow_size=2)

h = q.gf_hamilton_flowed_node(gf0, ga, fi)
print(f"Flowed Hamiltonian (node-local) = {h}")

q.end_with_mpi()
```

### Display Flow Info

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

fi = q.FlowInfo()
fi.add_flow(0, 0, 0.1, flow_size=3)
fi.add_flow(0, 1, 0.05, flow_size=2)
info_str = fi.show()
print(info_str)

q.end_with_mpi()
```
