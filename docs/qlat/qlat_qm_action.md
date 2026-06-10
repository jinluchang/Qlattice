# `qlat.qm_action` — Quantum-Mechanical Action for HMC

Source: `qlat/qlat/qm_action.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [The `QMAction` Class](#the-qmaction-class)
   - [Constructor](#constructor)
   - [Parameter Accessors](#parameter-accessors)
   - [Potential and Force](#potential-and-force)
   - [HMC Support](#hmc-support)
3. [Examples](#examples)

---

## Overview

The `qlat.qm_action` module provides `QMAction`, a Python wrapper around the
C-level quantum-mechanical action used in Hamiltonian Monte Carlo (HMC)
simulations. The action defines a confining potential `V(x, t)` with
parameters for barrier strength, finite-volume (FV) offsets, and measurement
windows.

The C implementation is created via `c.mk_qm_action` and freed on object
deletion via `c.free_qm_action`. All computation methods (potential, force,
field evolution) delegate to the corresponding `c.*` functions.

This module is typically used in the HMC update loop to:

1. Compute the action for Metropolis accept/reject (`action_node`).
2. Set conjugate momenta (`hmc_set_rand_momentum`).
3. Evolve fields and momenta (`hmc_field_evolve`, `hmc_set_force`).
4. Compute the Hamiltonian contribution from momenta (`hmc_m_hamilton_node`).

---

## The `QMAction` Class

### Constructor

```python
QMAction(
    alpha,
    beta,
    V_FV_min,
    FV_offset,
    TV_offset,
    barrier_strength,
    L,
    M,
    epsilon,
    t_FV_out,
    t_FV_mid,
    dt,
    measure_offset_L,
    measure_offset_M,
)
```

Create a new `QMAction` instance. All parameters are passed directly to the
C-level `mk_qm_action`.

| Parameter | Type | Description |
|---|---|---|
| `alpha` | `float` | Coupling constant |
| `beta` | `float` | Inverse temperature / coupling |
| `V_FV_min` | `float` | Minimum potential in the FV region |
| `FV_offset` | `float` | Finite-volume offset |
| `TV_offset` | `float` | Temporal-volume offset |
| `barrier_strength` | `float` | Strength of the confining barrier |
| `L` | `float` | Spatial extent parameter |
| `M` | `float` | Mass parameter |
| `epsilon` | `float` | Step-size parameter |
| `t_FV_out` | `int` | Outer FV time extent |
| `t_FV_mid` | `int` | Middle FV time extent |
| `dt` | `float` | Time step |
| `measure_offset_L` | `bool` | Measurement offset (L) |
| `measure_offset_M` | `bool` | Measurement offset (M) |

---

### Parameter Accessors

Each accessor returns the corresponding parameter stored in the C object.

| Method | Returns | Description |
|---|---|---|
| `alpha()` | `float` | Coupling constant |
| `beta()` | `float` | Inverse temperature / coupling |
| `barrier_strength()` | `float` | Confining barrier strength |
| `M()` | `float` | Mass parameter |
| `L()` | `float` | Spatial extent parameter |
| `t_FV_out()` | `int` | Outer FV time extent |
| `t_FV_mid()` | `int` | Middle FV time extent |
| `t_FV()` | `int` | Total FV time: `2 * t_FV_out + t_FV_mid` |
| `dt()` | `float` | Time step |

---

### Potential and Force

#### `V`

```python
V(x, t) -> float
```

Evaluate the confining potential at position `x = (x0, x1)` and time `t`.

| Parameter | Type | Description |
|---|---|---|
| `x` | array-like | Position; `x[0]` and `x[1]` are used |
| `t` | `int` | Time coordinate |
| **Returns** | `float` | Potential value |

---

#### `dV`

```python
dV(x, t) -> float
```

Evaluate the derivative of the potential with respect to `x[0]` at position
`x = (x0, x1)` and time `t`. (The underlying C function supports an
additional `idx` parameter to select the component, but the Python wrapper
always uses `idx=0`.)

| Parameter | Type | Description |
|---|---|---|
| `x` | array-like | Position; `x[0]` and `x[1]` are used |
| `t` | `int` | Time coordinate |
| **Returns** | `float` | `dV/dx0` |

---

### HMC Support

#### `action_node`

```python
action_node(f: FieldBase) -> float
```

Compute the local (per-node) action contribution for field `f`. Used in the
Metropolis accept/reject step.

---

#### `hmc_m_hamilton_node`

```python
hmc_m_hamilton_node(m: FieldBase) -> float
```

Compute the local kinetic-energy contribution from conjugate momenta `m`.

---

#### `sum_sq`

```python
sum_sq(f: FieldBase) -> float
```

Compute the sum of squares of field `f` (local node contribution).

---

#### `hmc_set_force`

```python
hmc_set_force(force: FieldBase, f: FieldBase) -> None
```

Compute the force (negative gradient of the action) for field `f` and store
it in `force`.

---

#### `hmc_field_evolve`

```python
hmc_field_evolve(f: FieldBase, m: FieldBase, step_size: float) -> None
```

Evolve field `f` by one leapfrog step: `f += step_size * m`. Modifies `f`
in place.

---

#### `hmc_set_rand_momentum`

```python
hmc_set_rand_momentum(m: FieldBase, rs: RngState) -> None
```

Set conjugate momenta `m` to random Gaussian values using the random number
state `rs`. Used at the beginning of each HMC trajectory.

---

## Examples

### Basic QMAction Setup

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

action = q.QMAction(
    alpha=1.0,
    beta=1.0,
    V_FV_min=0.0,
    FV_offset=0.0,
    TV_offset=0.0,
    barrier_strength=1.0,
    L=4.0,
    M=1.0,
    epsilon=0.1,
    t_FV_out=1,
    t_FV_mid=2,
    dt=0.01,
    measure_offset_L=False,
    measure_offset_M=False,
)

print(f"alpha = {action.alpha()}")
print(f"t_FV  = {action.t_FV()}")
print(f"V(0,0,0) = {action.V([0.0, 0.0], 0)}")

q.end_with_mpi()
```

### Evaluating Potential and Gradient

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

action = q.QMAction(
    alpha=1.0, beta=1.0, V_FV_min=0.0, FV_offset=0.0, TV_offset=0.0,
    barrier_strength=1.0, L=4.0, M=1.0, epsilon=0.1,
    t_FV_out=1, t_FV_mid=2, dt=0.01,
    measure_offset_L=False, measure_offset_M=False,
)

x = [1.0, 2.0]
t = 1
v = action.V(x, t)
dv = action.dV(x, t)
print(f"V({x}, {t}) = {v}")
print(f"dV({x}, {t}) = {dv}")

q.end_with_mpi()
```

### Copying a QMAction

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

action1 = q.QMAction(
    alpha=1.0, beta=1.0, V_FV_min=0.0, FV_offset=0.0, TV_offset=0.0,
    barrier_strength=1.0, L=4.0, M=1.0, epsilon=0.1,
    t_FV_out=1, t_FV_mid=2, dt=0.01,
    measure_offset_L=False, measure_offset_M=False,
)

# Copy via @= (target must be a valid QMAction instance)
action2 = q.QMAction(
    alpha=0.0, beta=0.0, V_FV_min=0.0, FV_offset=0.0, TV_offset=0.0,
    barrier_strength=0.0, L=0.0, M=0.0, epsilon=0.0,
    t_FV_out=0, t_FV_mid=0, dt=1.0,
    measure_offset_L=False, measure_offset_M=False,
)
action2 @= action1
print(f"Copied alpha = {action2.alpha()}")

q.end_with_mpi()
```
