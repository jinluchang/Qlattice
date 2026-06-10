# `qlat.hmc_stats` — HMC Diagnostic and Analysis Utilities

Source: `qlat/qlat/hmc_stats.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Functions](#functions)
   - [Force Magnitude Inspection](#force-magnitude-inspection)
   - [Gauge-Field Info with Wilson Flow](#gauge-field-info-with-wilson-flow)
3. [Examples](#examples)

---

## Overview

`hmc_stats` provides diagnostic utilities for monitoring and analysing
Hybrid Monte Carlo (HMC) trajectories. The functions wrap C/C++ routines
exposed through `qlat.c` and are typically called during or after HMC
evolution to:

- Inspect the magnitude distribution of the molecular-dynamics force
  across lattice sites.
- Save accumulated force-magnitude data to files for post-processing.
- Generate gauge-field information tables augmented with Wilson-flow
  energy density measurements.

All functions depend on the HMC machinery from `qlat.hmc` (imported via
`from qlat.hmc import *`).

---

## Functions

### Force Magnitude Inspection

#### `get_gm_force_magnitudes(gm_force, n_elems)`

```python
get_gm_force_magnitudes(gm_force: GaugeMomentum, n_elems: int) -> list
```

Compute and return the magnitude distribution of the gauge-momentum force
field. The force magnitude is binned into `n_elems` categories.

| Parameter | Type | Description |
|---|---|---|
| `gm_force` | `GaugeMomentum` | The molecular-dynamics force field |
| `n_elems` | `int` | Number of magnitude bins |

Returns a list of magnitude values.

#### `display_gm_force_magnitudes(gm_force, n_elems)`

```python
display_gm_force_magnitudes(gm_force: GaugeMomentum, n_elems: int) -> None
```

Compute the force magnitude distribution and display it to the log output.
Internally calls the C function and prints the result.

| Parameter | Type | Description |
|---|---|---|
| `gm_force` | `GaugeMomentum` | The molecular-dynamics force field |
| `n_elems` | `int` | Number of magnitude bins |

#### `save_gm_force_magnitudes_list(fn)`

```python
save_gm_force_magnitudes_list(fn: str) -> None
```

Save the accumulated list of gauge-momentum force magnitudes to a file.
Creates parent directories if needed.

| Parameter | Type | Description |
|---|---|---|
| `fn` | `str` | Output file path |

---

### Gauge-Field Info with Wilson Flow

#### `display_gauge_field_info_table_with_wilson_flow(fn_gf_info, fn_wilson_flow_energy, gf, flow_time, flow_steps, steps, c1=0.0)`

```python
display_gauge_field_info_table_with_wilson_flow(
    fn_gf_info: str,
    fn_wilson_flow_energy: str,
    gf: GaugeField,
    flow_time: float,
    flow_steps: int,
    steps: int,
    c1: float = 0.0,
) -> None
```

Generate a gauge-field information table with Wilson-flow energy density
measurements. Applies Wilson flow to `gf` and records both plaquette-based
observables and the flowed energy density at each step.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fn_gf_info` | `str` | — | Output path for gauge-field info table |
| `fn_wilson_flow_energy` | `str` | — | Output path for Wilson-flow energy data |
| `gf` | `GaugeField` | — | Gauge field to analyse |
| `flow_time` | `float` | — | Total Wilson flow time |
| `flow_steps` | `int` | — | Number of flow steps |
| `steps` | `int` | — | Number of measurement intervals |
| `c1` | `float` | `0.0` | Symanzik improvement coefficient |

---

## Examples

### Inspecting Force Magnitudes During HMC

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
rs = q.RngState("force-test")

gf = q.GaugeField(geo)
gf.set_unit()

ga = q.GaugeAction(6.0)
gm = q.GaugeMomentum(geo)
q.set_rand_gauge_momentum(gm, 1.0, rs.split("gm"))

q.set_gm_force(gm, gf, ga)
q.display_gm_force_magnitudes(gm, 5)

q.end_with_mpi()
```

### Wilson-Flow Analysis After HMC Trajectories

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
rs = q.RngState("hmc-wflow")

gf = q.GaugeField(geo)
gf.set_unit()

ga = q.GaugeAction(5.5)

for traj in range(3):
    q.run_hmc_pure_gauge(
        gf, ga, traj, rs.split(f"traj-{traj}"), is_always_accept=True
    )
    plaq = q.gf_avg_plaq(gf)
    print(f"traj={traj}  plaq={plaq:.12E}")

    q.display_gauge_field_info_table_with_wilson_flow(
        f"results/gf_info/traj={traj}.lat",
        f"results/wilson_flow_energy/traj={traj}.lat",
        gf,
        0.1,
        5,
        2,
    )

q.end_with_mpi()
```

### Saving Force Magnitude Data

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
rs = q.RngState("save-force")

gf = q.GaugeField(geo)
gf.set_unit()

ga = q.GaugeAction(6.0)
gm = q.GaugeMomentum(geo)

for traj in range(3):
    q.run_hmc_pure_gauge(
        gf, ga, traj, rs.split(f"traj-{traj}"), is_always_accept=True
    )
    q.set_gm_force(gm, gf, ga)
    q.display_gm_force_magnitudes(gm, 5)

q.save_gm_force_magnitudes_list("results/gm_force_info/all.dat")

q.end_with_mpi()
```
