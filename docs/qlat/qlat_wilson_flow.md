# `qlat.wilson_flow` — Wilson Flow and Stout Smearing Utilities

Source: `qlat/qlat/wilson_flow.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Energy Density Functions](#energy-density-functions)
   - [`gf_energy_density`](#gf_energy_density)
   - [`gf_energy_density_field`](#gf_energy_density_field)
   - [`gf_energy_density_dir_field`](#gf_energy_density_dir_field)
3. [Wilson Flow](#wilson-flow)
   - [`gf_wilson_flow_force`](#gf_wilson_flow_force)
   - [`gf_wilson_flow_step`](#gf_wilson_flow_step)
   - [`gf_wilson_flow`](#gf_wilson_flow)
   - [`gf_energy_derivative_density_field`](#gf_energy_derivative_density_field)
   - [`gf_plaq_flow_force`](#gf_plaq_flow_force)
4. [Stout Smearing](#stout-smearing)
   - [`gf_stout_smear`](#gf_stout_smear)
   - [`gf_block_stout_smear`](#gf_block_stout_smear)
   - [`gf_local_stout_smear`](#gf_local_stout_smear)
   - [`gf_local_avg_plaq`](#gf_local_avg_plaq)
5. [Local Tree Gauge Fixing](#local-tree-gauge-fixing)
   - [`mk_local_tree_gauge_f_dir`](#mk_local_tree_gauge_f_dir)
   - [`gt_local_tree_gauge`](#gt_local_tree_gauge)
   - [`gt_block_tree_gauge`](#gt_block_tree_gauge)
6. [Examples](#examples)

---

## Overview

`wilson_flow` implements the Yang-Mills gradient flow (Wilson flow) and
related gauge field operations used in lattice QCD simulations. The Wilson
flow evolves a gauge field in a fictitious flow time `t`, smoothing
short-distance fluctuations while preserving long-distance topology. This is
essential for computing the running coupling, energy density, and topological
charge.

The module also provides stout smearing and a local tree gauge fixing
procedure that can be used for gauge fixing within spatial blocks.

Key references:
- [arXiv:1006.4518](https://arxiv.org/abs/1006.4518) — Wilson flow and energy density
- [arXiv:1203.4469](https://arxiv.org/abs/1203.4469) — Topological susceptibility from Wilson flow

---

## Energy Density Functions

### `gf_energy_density`

```python
gf_energy_density(gf: GaugeField) -> float
```

Compute the total plaquette energy density (Euclidean) of the gauge field.
Returns the global sum of `E(x) = -1/2 * Tr(F_{mu,nu} F_{mu,nu})` over all
sites. The result is not normalized by volume.

### `gf_energy_density_field`

```python
gf_energy_density_field(gf: GaugeField) -> FieldRealD
```

Return a `FieldRealD` with `multiplicity == 1` containing the local energy
density at each lattice site. This is the per-site version of
`gf_energy_density`.

### `gf_energy_density_dir_field`

```python
gf_energy_density_dir_field(gf: GaugeField) -> FieldRealD
```

Return a `FieldRealD` with `multiplicity == 6` containing the energy density
decomposed by the six Lorentz directions (`(mu, nu)` with `mu < nu`).

For a smeared gauge field the relation to the plaquette is approximately:

```
gf_energy_density_dir_field(gf)[:]  ~  6 * (1 - gf_plaq_field(gf)[:])
```

Summing over the last axis recovers `gf_energy_density_field`:

```
gf_energy_density_dir_field(gf)[:].sum(-1, keepdims=True)
    == gf_energy_density_field(gf)[:]
```

See [arXiv:1006.4518 Eq. (2.1)](https://arxiv.org/abs/1006.4518) and
[arXiv:1203.4469](https://arxiv.org/abs/1203.4469).

---

## Wilson Flow

### `gf_wilson_flow_force`

```python
gf_wilson_flow_force(gf: GaugeField, c1: float = 0.0) -> GaugeMomentum
```

Compute the Wilson flow force (the Lie algebra valued derivative of the
action) for gauge field `gf`. Parameter `c1` controls the Symanzik
improvement (`c1 = 0` for the standard Wilson plaquette action). Returns
the force as a `GaugeMomentum`.

### `gf_wilson_flow_step`

```python
gf_wilson_flow_step(
    gf: GaugeField,
    epsilon: float,
    *,
    c1: float = 0.0,
    wilson_flow_integrator_type: str = None,
)
```

Evolve the gauge field `gf` in place by one Wilson flow step of size
`epsilon`. The default integrator is the 3rd-order Runge-Kutta scheme from
[arXiv:1006.4518](https://arxiv.org/abs/1006.4518). An Euler integrator is
also available by passing `wilson_flow_integrator_type="euler"`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Gauge field (modified in place) |
| `epsilon` | `float` | — | Flow step size |
| `c1` | `float` | `0.0` | Symanzik improvement parameter |
| `wilson_flow_integrator_type` | `str` or `None` | `None` | `"runge-kutta"` (default) or `"euler"` |

### `gf_wilson_flow`

```python
gf_wilson_flow(
    gf: GaugeField,
    flow_time: float,
    steps: int,
    *,
    c1: float = 0.0,
    existing_flow_time: float = 0.0,
    wilson_flow_integrator_type: str = None,
) -> list[float]
```

Evolve the gauge field `gf` in place by a total Wilson flow time
`flow_time`, divided into `steps` equal sub-steps. Returns a list of energy
density values recorded after each sub-step. The printed `t^2 E` values are
used to determine the scale `t_0` (defined by `t^2 E = 0.3`).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Gauge field (modified in place) |
| `flow_time` | `float` | — | Total flow time |
| `steps` | `int` | — | Number of sub-steps |
| `c1` | `float` | `0.0` | Symanzik improvement parameter |
| `existing_flow_time` | `float` | `0.0` | Offset for reporting (e.g., if field was already flowed) |
| `wilson_flow_integrator_type` | `str` or `None` | `None` | `"runge-kutta"` or `"euler"` |

### `gf_energy_derivative_density_field`

```python
gf_energy_derivative_density_field(
    gf: GaugeField,
    *,
    epsilon: float = 0.0125,
    c1: float = 0.0,
    wilson_flow_integrator_type: str = None,
) -> FieldRealD
```

Compute the derivative `dE/dt` of the energy density with respect to flow
time at `t = 0` using a symmetric finite difference of step `2 * epsilon`.
Returns a `FieldRealD` with `multiplicity == 1`. This is related to the
topological charge via the Yang-Mills gradient flow renormalization.

### `gf_plaq_flow_force`

```python
gf_plaq_flow_force(
    gf: GaugeField,
    plaq_factor: FieldRealD,
) -> GaugeMomentum
```

Compute the flow force with a plaquette-dependent beta factor, allowing
spatially varying couplings. `plaq_factor` must have `multiplicity == 6`
(one factor per plaquette direction). The ordering of plaquettes matches
`gf_plaq_field`.

---

## Stout Smearing

### `gf_stout_smear`

```python
gf_stout_smear(
    gf: GaugeField,
    step_size: float,
    num_step: int = 1,
    *,
    method: str = None,
)
```

Apply stout smearing to `gf` in place for `num_step` steps. `method`
selects the implementation: `None` or `"force"` calls
`gf_wilson_flow_force` + `gf_evolve` directly (default), `"stout"` uses
`gf_block_stout_smear`, `"wilson-flow"` uses `gf_wilson_flow_step` with
the Euler integrator and `c1=0`. All three are mathematically equivalent.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Gauge field (modified in place) |
| `step_size` | `float` | — | Smearing strength ρ (typical: 0.1–0.125) |
| `num_step` | `int` | `1` | Number of smearing iterations |
| `method` | `str` | `None` | `None`/`"stout"`, `"wilson-flow"`, or `"force"` |

### `gf_block_stout_smear`

```python
gf_block_stout_smear(
    gf: GaugeField,
    block_site: Coordinate,
    step_size: float,
)
```

Apply stout smearing to `gf` in place. If `block_site` is an empty
`Coordinate()`, smearing is applied globally. Otherwise, each spatial block
of size `block_site` is smeared independently (useful for local gauge
fixing workflows).

### `gf_local_stout_smear`

```python
gf_local_stout_smear(
    gf: GaugeField,
    block_site: Coordinate,
    step_size: float,
)
```

Apply stout smearing locally (no MPI communication). If `block_site` is an
empty `Coordinate()`, no blocking is applied. Otherwise, each block is
smeared independently.

### `gf_local_avg_plaq`

```python
gf_local_avg_plaq(
    gf: GaugeField,
    block_site: Coordinate,
) -> float
```

Compute the local (no global sum) average plaquette within blocks of size
`block_site`. Used as a diagnostic in the `gt_local_tree_gauge` workflow.
If `block_site == Coordinate()`, defaults to `geo.node_site`.

---

## Local Tree Gauge Fixing

### `mk_local_tree_gauge_f_dir`

```python
mk_local_tree_gauge_f_dir(
    geo: Geometry,
    block_site: Coordinate,
    is_uniform: bool,
    rs: RngState,
) -> FieldInt
```

Generate a `FieldInt` (multiplicity 1) specifying the tree direction at each
lattice site. If `f_dir[x] == 5` for a site, that site is a tree root (the
gauge transformation is the identity). Used as input to
`gt_local_tree_gauge`.

### `gt_local_tree_gauge`

```python
gt_local_tree_gauge(
    gf: GaugeField,
    f_dir: FieldInt,
    num_step: int,
) -> GaugeTransform
```

Compute a gauge transformation `gt_inv` using a local tree gauge fixing
procedure. The fixed gauge field is obtained as:

```python
gf_fixed = gt_inv.inv() * gf
```

Note: this convention is the **inverse** of the usual `gt * gf` convention
used elsewhere.

### `gt_block_tree_gauge`

```python
gt_block_tree_gauge(
    gf: GaugeField,
    *,
    block_site: Coordinate = None,
    new_size_node: Coordinate = None,
    is_uniform: bool = True,
    stout_smear_step_size: float = 0.125,
    num_smear_step: int = 4,
    f_dir_list: list = None,
    rs_f_dir: RngState = None,
) -> tuple[GaugeTransform, list[FieldInt]]
```

High-level block tree gauge fixing. Returns `(gt_inv, f_dir_list)`.

Steps performed:
1. Optionally reshuffle the field to match `new_size_node`.
2. Apply `num_smear_step` rounds of local stout smearing to each block.
3. Compute local tree gauge fixing for each block.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field (not modified) |
| `block_site` | `Coordinate` | `[4,4,4,4]` | Block dimensions |
| `new_size_node` | `Coordinate` | Derived from `total_site / block_site` | MPI node layout |
| `is_uniform` | `bool` | `True` | Whether tree directions are uniform |
| `stout_smear_step_size` | `float` | `0.125` | Stout smearing step size |
| `num_smear_step` | `int` | `4` | Number of smearing steps |
| `f_dir_list` | `list` or `None` | `None` | Pre-computed tree directions; if `None`, generated from `rs_f_dir` |
| `rs_f_dir` | `RngState` or `None` | `None` | RNG state for generating `f_dir_list` |

---

## Examples

### Wilson Flow and Scale Setting

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
rng = q.RngState("seed")
gf.set_rand(rng)

flow_time = 1.0
steps = 100
energy_list = q.gf_wilson_flow(gf, flow_time, steps)

for i, e in enumerate(energy_list):
    t = (i + 1) * (flow_time / steps)
    print(f"t={t:.4f}  E={e:.6f}  t^2*E={t*t*e:.6f}")

q.end_with_mpi()
```

### Single Wilson Flow Step

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
rng = q.RngState("seed")
gf.set_rand(rng)

e_before = q.gf_energy_density(gf)
q.gf_wilson_flow_step(gf, epsilon=0.01)
e_after = q.gf_energy_density(gf)

print(f"E before flow: {e_before}")
print(f"E after flow:  {e_after}")

q.end_with_mpi()
```

### Energy Density Direction Field

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
rng = q.RngState("seed")
gf.set_rand(rng)

ed_dir = q.gf_energy_density_dir_field(gf)
ed_total = q.gf_energy_density_field(gf)

import numpy as np
np.testing.assert_allclose(
    np.asarray(ed_dir).sum(axis=-1, keepdims=True),
    np.asarray(ed_total),
    atol=1e-12,
)
print("Direction decomposition matches total energy density.")

q.end_with_mpi()
```

### Block Stout Smearing

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([8, 8, 8, 8])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
rng = q.RngState("seed")
gf.set_rand(rng)

block_site = q.Coordinate([4, 4, 4, 4])
step_size = 0.1

print(f"Plaq before smear: {q.gf_avg_plaq(gf)}")
q.gf_block_stout_smear(gf, block_site, step_size)
print(f"Plaq after smear:  {q.gf_avg_plaq(gf)}")

q.end_with_mpi()
```

### Block Tree Gauge Fixing

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([8, 8, 8, 8])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
rng = q.RngState("seed")
gf.set_rand(rng)

block_site = q.Coordinate([4, 4, 4, 4])
rs = q.RngState("seed-gt")

gt_inv, f_dir_list = q.gt_block_tree_gauge(
    gf,
    block_site=block_site,
    is_uniform=True,
    rs_f_dir=rs,
)

gf_fixed = gt_inv.inv() * gf
print(f"Plaq before fix: {q.gf_avg_plaq(gf)}")
print(f"Plaq after fix:  {q.gf_avg_plaq(gf_fixed)}")

q.end_with_mpi()
```
