# `qlat.gauge_action` — Gauge Action Parameters

Source: `qlat/qlat/gauge_action.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`GaugeAction` Class](#gaugeaction-class)
   - [Constructor](#constructor)
   - [Methods](#methods)
3. [Examples](#examples)

---

## Overview

`gauge_action` defines the `GaugeAction` class, a lightweight container for the
parameters that specify a lattice gauge action. In the standard Wilson
formulation the action is parameterised by a single coupling `beta`; adding a
rectangle term introduces a second parameter `c1` (improved actions such as
Iwasaki or DBW2).

The stored parameters are consumed by the HMC routines in `qlat.hmc` (e.g.
`gf_hamilton_node`, `set_gm_force`) to evaluate the gauge Hamiltonian and its
molecular-dynamics force.

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

ga = q.GaugeAction(5.5)          # Wilson action, beta=5.5
print(ga.beta())                 # 5.5
print(ga.c1())                   # 0.0

q.end_with_mpi()
```

---

## `GaugeAction` Class

### Constructor

```python
GaugeAction(beta: float, c1: float = 0.0)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `beta` | `float` | — | Inverse gauge coupling (6 / g^2) |
| `c1` | `float` | `0.0` | Coefficient of the rectangle (improvement) term |

- `c1 = 0.0` gives the standard Wilson plaquette action.
- `c1 = -0.331` corresponds to the Iwasaki action (the C++ default).

```python
ga_wilson  = q.GaugeAction(6.0)            # pure Wilson
ga_iwasaki = q.GaugeAction(2.62, -0.331)   # Iwasaki improved
```

---

### Methods

#### `beta() -> float`

Return the inverse gauge coupling parameter.

#### `c1() -> float`

Return the rectangle-term coefficient.

---

## Examples

### Basic Usage

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

ga = q.GaugeAction(5.5, -0.331)
print(f"beta = {ga.beta()}")     # 5.5
print(f"c1   = {ga.c1()}")       # -0.331

# Copy via @=
ga2 = q.GaugeAction(0.0)
ga2 @= ga
print(f"ga2 beta = {ga2.beta()}")  # 5.5

q.end_with_mpi()
```

### Gauge Hamiltonian with `GaugeAction`

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
gf.set_rand(rs.split("gf-init"), 0.5, 10)
gf.unitarize()

ga = q.GaugeAction(6.0)
h = q.gf_hamilton_node(gf, ga)
print(f"gauge Hamiltonian (node) = {h}")

q.end_with_mpi()
```
