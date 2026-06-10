# `qlat.hmc` — Hybrid Monte Carlo for Pure-Gauge Simulations

Source: `qlat/qlat/hmc.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`GaugeMomentum` Class](#gaugemomentum-class)
   - [Constructor](#constructor)
   - [Methods](#methods)
3. [Free Functions](#free-functions)
   - [Momentum Initialization](#momentum-initialization)
   - [Hamiltonian Helpers](#hamiltonian-helpers)
   - [Molecular-Dynamics Evolution](#molecular-dynamics-evolution)
   - [Force Computation](#force-computation)
   - [Gauge-Transformation Utilities](#gauge-transformation-utilities)
   - [Matrix Basis Conversions](#matrix-basis-conversions)
   - [Metropolis Accept/Reject](#metropolis-acceptreject)
   - [HMC Drivers](#hmc-drivers)
4. [Examples](#examples)

---

## Overview

`hmc` provides the building blocks for a Hybrid Monte Carlo (HMC) simulation of
pure SU(3) lattice gauge theory. The key components are:

- **`GaugeMomentum`** — a `FieldColorMatrix` with multiplicity 4 (one
  anti-Hermitian matrix per spacetime direction) representing the conjugate
  momentum to the gauge field.
- **Hamiltonian evaluation** — `gm_hamilton_node` (kinetic energy) and
  `gf_hamilton_node` (potential energy, using a `GaugeAction`).
- **Molecular-dynamics evolution** — `gf_evolve` / `gf_evolve_dual` advance the
  gauge field by a leapfrog step; `gf_evolve_fa` / `gf_evolve_fa_dual` support
  field-dependent (flavour-action) weights.
- **Force computation** — `set_gm_force` computes the derivative of the gauge
  Hamiltonian with respect to the gauge field.
- **Metropolis step** — `metropolis_accept` decides whether to accept the
  proposed configuration based on the Hamiltonian change.
- **High-level driver** — `run_hmc_pure_gauge` performs a complete HMC
  trajectory: generate momenta, integrate, accept/reject.

All functions operate in a distributed MPI environment; `_node` suffixes
indicate node-local (non-communicating) quantities.

---

## `GaugeMomentum` Class

`GaugeMomentum` is a subclass of `FieldColorMatrix` with fixed multiplicity 4
(one anti-Hermitian 3x3 matrix per direction). It serves as the conjugate
momentum field in HMC.

### Constructor

```python
GaugeMomentum(geo: Geometry = None)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `geo` | `Geometry` | `None` | Lattice geometry; if `None`, the field is uninitialized |

### Methods

#### `set_rand(rng: RngState, sigma: float = 1.0)`

Fill the momentum field with Gaussian random anti-Hermitian matrices drawn
from `rng` with standard deviation `sigma`.

#### `set_rand_fa(mf: FieldRealD, rng: RngState)`

Fill the momentum field with Gaussian random anti-Hermitian matrices using a
per-site, per-direction weight from `mf` (multiplicity 4). This is used for
flavour-action (field-dependent) HMC.

---

## Free Functions

### Momentum Initialization

#### `set_rand_gauge_momentum(gm, sigma, rng)`

```python
set_rand_gauge_momentum(gm: GaugeMomentum, sigma: float, rng: RngState) -> None
```

Fill `gm` with Gaussian random anti-Hermitian matrices (standard deviation
`sigma`).

#### `set_rand_gauge_momentum_fa(gm, mf, rng)`

```python
set_rand_gauge_momentum_fa(gm: GaugeMomentum, mf: FieldRealD, rng: RngState) -> None
```

Flavour-action variant: weight each direction by the corresponding element of
`mf`.

---

### Hamiltonian Helpers

#### `gm_hamilton_node(gm) -> float`

```python
gm_hamilton_node(gm: GaugeMomentum) -> float
```

Return the node-local kinetic energy: sum of `Tr[p^2] / 2` over all sites and
directions.

#### `gm_hamilton_node_fa(gm, mf) -> float`

```python
gm_hamilton_node_fa(gm: GaugeMomentum, mf: FieldRealD) -> float
```

Flavour-action variant of the kinetic energy, weighted by `mf`.

#### `gf_hamilton_node(gf, ga) -> float`

```python
gf_hamilton_node(gf: GaugeField, ga: GaugeAction) -> float
```

Return the node-local gauge Hamiltonian (plaquette and optional rectangle
terms) using the parameters in `ga`.

---

### Molecular-Dynamics Evolution

#### `gf_evolve(gf, gm, step_size)`

```python
gf_evolve(gf: GaugeField, gm: GaugeMomentum, step_size: float) -> None
```

Update the gauge field: `U <- exp(i * step_size * P) * U` for each link.

#### `gf_evolve_dual(gf, gm_dual, step_size)`

```python
gf_evolve_dual(gf: GaugeField, gm_dual: GaugeMomentum, step_size: float) -> None
```

Dual (adjoint) evolution: `U <- U * exp(-i * step_size * P_dual)`.

#### `gf_evolve_fa(gf, gm, mf, step_size)`

```python
gf_evolve_fa(gf: GaugeField, gm: GaugeMomentum, mf: FieldRealD, step_size: float) -> None
```

Flavour-action variant of `gf_evolve`, weighted per-site/direction by `mf`.

#### `gf_evolve_fa_dual(gf, gm_dual, mf_dual, step_size)`

```python
gf_evolve_fa_dual(gf: GaugeField, gm_dual: GaugeMomentum, mf_dual: FieldRealD, step_size: float) -> None
```

Flavour-action dual evolution.

---

### Force Computation

#### `set_gm_force(gm_force, gf, ga)`

```python
set_gm_force(gm_force: GaugeMomentum, gf: GaugeField, ga: GaugeAction) -> None
```

Compute the molecular-dynamics force on the gauge field and store the result
in `gm_force`.

#### `set_gm_force_dual(gm_force_dual, gf, gm_force)`

```python
set_gm_force_dual(gm_force_dual: GaugeMomentum, gf: GaugeField, gm_force: GaugeMomentum) -> None
```

Compute the dual force contribution.

---

### Gauge-Transformation Utilities

#### `project_gauge_transform(gm, gm_dual, mf, mf_dual) -> float`

```python
project_gauge_transform(gm: GaugeMomentum, gm_dual: GaugeMomentum,
                        mf: FieldRealD, mf_dual: FieldRealD) -> float
```

Project out the pure gauge-transformation component from the momentum fields.
Returns the projected norm.

#### `set_gauge_transform_momentum(gm, gm_dual, gtm)`

```python
set_gauge_transform_momentum(gm: GaugeMomentum, gm_dual: GaugeMomentum,
                             gtm: FieldColorMatrix) -> None
```

Set `gm` and `gm_dual` so that a combined evolution with `gf_evolve` and
`gf_evolve_dual` produces a pure gauge transformation parameterised by `gtm`.

---

### Matrix Basis Conversions

#### `set_anti_hermitian_matrix_from_basis(fc, basis)`

```python
set_anti_hermitian_matrix_from_basis(fc: FieldColorMatrix, basis: FieldRealD) -> None
```

Convert a real-valued basis representation to anti-Hermitian matrices:

```
fc[x, m] = sum_a  T_a * basis[x, m * 8 + a]
```

where `T_a` are the 8 SU(3) generators satisfying `T_a^dag = -T_a` and
`Tr[T_a T_b] = -2 delta_ab`.

#### `set_basis_from_anti_hermitian_matrix(basis, fc)`

```python
set_basis_from_anti_hermitian_matrix(basis: FieldRealD, fc: FieldColorMatrix) -> None
```

Inverse of `set_anti_hermitian_matrix_from_basis`: extract the 8-component
real basis from each anti-Hermitian matrix.

---

### Metropolis Accept/Reject

#### `metropolis_accept(delta_h, traj, rs) -> (bool, float)`

```python
metropolis_accept(delta_h: float, traj: int, rs: RngState) -> (bool, float)
```

Perform the Metropolis accept/reject step. The proposed configuration is
accepted with probability `min(1, exp(-delta_h))`.

| Parameter | Type | Description |
|---|---|---|
| `delta_h` | `float` | Change in Hamiltonian (H_new - H_old) |
| `traj` | `int` | Trajectory number (for logging) |
| `rs` | `RngState` | Random state for the accept/reject draw |

Returns `(flag, accept_prob)` where `flag` is `True` if accepted and
`accept_prob` is the acceptance probability.

---

### HMC Drivers

#### `gm_evolve_fg_pure_gauge(gm, gf_init, ga, fg_dt, dt)`

```python
gm_evolve_fg_pure_gauge(gm, gf_init, ga, fg_dt, dt) -> None
```

First-order (force-gradient) momentum update for pure-gauge HMC. Used
internally by `run_hmc_evolve_pure_gauge`.

| Parameter | Type | Description |
|---|---|---|
| `gm` | `GaugeMomentum` | Momentum field (modified in-place) |
| `gf_init` | `GaugeField` | Initial gauge field (not modified) |
| `ga` | `GaugeAction` | Gauge action parameters |
| `fg_dt` | `float` | Force-gradient step size |
| `dt` | `float` | Leapfrog step size |

#### `run_hmc_evolve_pure_gauge(gm, gf, ga, rs, n_step, md_time=1.0) -> float`

```python
run_hmc_evolve_pure_gauge(gm, gf, ga, rs, n_step, md_time=1.0) -> float
```

Run the molecular-dynamics trajectory for pure-gauge HMC using a
4th-order Omelyan integrator (with force-gradient improvement). Both `gm`
and `gf` are modified in-place.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gm` | `GaugeMomentum` | — | Conjugate momentum field |
| `gf` | `GaugeField` | — | Gauge field to evolve |
| `ga` | `GaugeAction` | — | Gauge action parameters |
| `rs` | `RngState` | — | Random state |
| `n_step` | `int` | — | Number of MD steps |
| `md_time` | `float` | `1.0` | Total molecular-dynamics time |

Returns the change in Hamiltonian `delta_h = H_final - H_initial`.

#### `run_hmc_pure_gauge(gf, ga, traj, rs, ...) -> (bool, float)`

```python
run_hmc_pure_gauge(gf, ga, traj, rs, *,
                   is_reverse_test=False, n_step=6, md_time=1.0,
                   is_always_accept=False) -> (bool, float)
```

Perform one complete HMC trajectory for pure-gauge theory:

1. Generate Gaussian momenta from `rs`.
2. Integrate the MD equations via `run_hmc_evolve_pure_gauge`.
3. (Optional) Verify time-reversibility.
4. Apply the Metropolis accept/reject step.
5. If accepted (or `is_always_accept`), update `gf` in-place.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Gauge field (updated if accepted) |
| `ga` | `GaugeAction` | — | Gauge action parameters |
| `traj` | `int` | — | Trajectory number |
| `rs` | `RngState` | — | Random state (split internally by `traj`) |
| `is_reverse_test` | `bool` | `False` | Run a reversibility check after the trajectory |
| `n_step` | `int` | `6` | Number of MD steps |
| `md_time` | `float` | `1.0` | Total MD time |
| `is_always_accept` | `bool` | `False` | Skip the Metropolis step and always accept |

Returns `(flag, delta_h)` where `flag` is `True` if the trajectory was
accepted and `delta_h` is the Hamiltonian change.

---

## Examples

### Basic HMC Trajectory

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

gf = q.GaugeField(geo)
gf.unitarize()

ga = q.GaugeAction(6.0)
rs = q.RngState("test-hmc")

flag, delta_h = q.run_hmc_pure_gauge(gf, ga, 0, rs, n_step=6, md_time=1.0)
print(f"accepted={flag}, delta_h={delta_h}")

q.end_with_mpi()
```

### Gauge Momentum Operations

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

gm = q.GaugeMomentum(geo)
rs = q.RngState("momentum-test")
gm.set_rand(rs, 1.0)

ke = q.gm_hamilton_node(gm)
print(f"kinetic energy (node) = {ke}")

q.end_with_mpi()
```

### Matrix Basis Conversion

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

fc = q.FieldColorMatrix(geo, 4)
basis = q.FieldRealD(geo, 4 * 8)

# Convert anti-Hermitian matrices to basis and back
q.set_basis_from_anti_hermitian_matrix(basis, fc)
q.set_anti_hermitian_matrix_from_basis(fc, basis)

q.end_with_mpi()
```

### Multiple Trajectories with Acceptance Monitoring

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

gf = q.GaugeField(geo)
gf.unitarize()

ga = q.GaugeAction(6.0)
rs = q.RngState("hmc-run")

n_traj = 10
accept_count = 0
for traj in range(n_traj):
    flag, delta_h = q.run_hmc_pure_gauge(gf, ga, traj, rs.split(f"{traj}"),
                                         n_step=6, md_time=1.0)
    if flag:
        accept_count += 1

print(f"acceptance rate = {accept_count / n_traj * 100:.1f}%")

q.end_with_mpi()
```
