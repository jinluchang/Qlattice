# `qlat.instanton_map` — Instanton Detection and Topological Charge Mapping

Source: `qlat/qlat/instanton_map.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Plaq-Dependent Flow](#plaq-dependent-flow)
   - [`gf_flow_topo`](#gf_flow_topo)
   - [`gf_energy_derivative_density_field_topo`](#gf_energy_derivative_density_field_topo)
3. [Plaquette Coordinate Utilities](#plaquette-coordinate-utilities)
   - [`mk_plaq_xg_arr`](#mk_plaq_xg_arr)
   - [`get_extreme_plaq_xg_list`](#get_extreme_plaq_xg_list)
   - [`get_group_extreme_plaq_xg_list`](#get_group_extreme_plaq_xg_list)
4. [Instanton List Processing](#instanton-list-processing)
   - [`process_inst_list`](#process_inst_list)
   - [`get_tot_topo_count`](#get_tot_topo_count)
   - [`get_tot_inst_count`](#get_tot_inst_count)
5. [Display Utilities](#display-utilities)
   - [`displayln_info_topo_info`](#displayln_info_topo_info)
   - [`displayln_info_p_inst`](#displayln_info_p_inst)
   - [`displayln_info_p_inst_list`](#displayln_info_p_inst_list)
6. [InstantonMap Class](#instantonmap-class)
   - [`InstantonMap.__init__`](#instantonmap__init__)
   - [`InstantonMap.shift`](#instantonmap_shift)
   - [`InstantonMap.convert_xg`](#instantonmap_convert_xg)
   - [`InstantonMap.half_lattice`](#instantonmap_half_lattice)
   - [`InstantonMap.acc_time`](#instantonmap_acc_time)
   - [`InstantonMap.acc_topo_info`](#instantonmap_acc_topo_info)
7. [Top-Level Computation](#top-level-computation)
   - [`compute_inst_map`](#compute_inst_map)
   - [`displayln_info_inst_map_obj`](#displayln_info_inst_map_obj)
   - [`smear_measure_topo`](#smear_measure_topo)
8. [Examples](#examples)

---

## Overview

`instanton_map` provides tools for detecting and tracking instantons (and
anti-instantons) in lattice gauge field configurations. The module uses
modified Wilson flow with plaquette-dependent couplings to either freeze,
shrink, or localize topological structures, then identifies instanton
candidates by grouping lattice sites with anomalously low plaquette values.

The key idea is that instantons manifest as localized regions where the
plaquette deviates significantly from unity. By flowing the gauge field
with specialized flow types, one can either preserve the topological
charge (freeze flow) or enhance tunnelling (shrink flow), enabling
robust instanton counting and topological charge estimation.

Standard Wilson flow action:
$$
S_\mathrm{Wilson} = \frac{\beta}{2}\sum_{x,\mu,\nu} \Big(1 - \frac{1}{3}\mathrm{Re}\mathrm{Tr} U_{\mu,\nu}\Big)
$$

Modified flow with plaquette-dependent coupling:
$$
S_f = -\frac{\beta}{2}\sum_{x,\mu,\nu} f\Big(\frac{1}{3}\mathrm{Re}\mathrm{Tr} U_{\mu,\nu}\Big)
$$

---

## Plaq-Dependent Flow

### `gf_flow_topo`

```python
gf_flow_topo(
    gf: GaugeField,
    step_size: float,
    flow_type: str | float | None = None,
    wilson_flow_integrator_type: str | None = None,
) -> float
```

Perform one step of plaquette-dependent flow on `gf` in place. Returns
`step_size / norm`, the effective step size used.

For standard Symanzik flows (`"Wilson"`, `"Iwasaki"`, `"DBW2"`, or a
numeric `c1` value), delegates to `gf_wilson_flow_step` with Euler
integrator by default.

For specialized flow types, a plaquette-dependent force is computed:

| Flow Type | Description | Parameters |
|---|---|---|
| `"Wilson"` | Standard Wilson flow (`c1 = 0.0`) | — |
| `"Iwasaki"` | Iwasaki action (`c1 = -0.331`) | — |
| `"DBW2"` | DBW2 action (`c1 = -1.4008`) | — |
| `"Freeze"` | Suppresses plaquette deviations; prevents topological tunnelling | — |
| `"Shrink"` | Shrinks instanton size; enhances topological tunnelling | `eps=0.005`, `b=0.5` |
| `"Localize"` | Shrinks large instantons, prevents small instanton tunnelling | `eps=0.002`, `b=50` |
| `"Preserve"` | Mimics Wilson flow, prevents small instanton tunnelling | `b=50` |

The freeze derivative: `df/dp = 1 - p`
The shrink derivative: `df/dp = eps / (1 - p + eps) + b`
The localize derivative: `df/dp = eps / (1 - p + eps) + b * (1 - p)`
The preserve derivative: `df/dp = max(1, b * (1 - p))`

### `gf_energy_derivative_density_field_topo`

```python
gf_energy_derivative_density_field_topo(
    gf: GaugeField,
    *,
    epsilon: float = 0.0125,
    flow_typw: str | float | None = None,
    wilson_flow_integrator_type: str | None = None,
) -> FieldRealD
```

Compute the derivative `dE/dt` of the energy density with respect to flow
time using a symmetric finite difference of step `2 * epsilon`. Supports
the same flow types as `gf_flow_topo`. Returns a `FieldRealD` with
multiplicity 1.

---

## Plaquette Coordinate Utilities

### `mk_plaq_xg_arr`

```python
mk_plaq_xg_arr(geo: Geometry) -> np.ndarray
```

Build an array of plaquette center coordinates for all local sites.
Returns an array of shape `(n_points, 6, 4)` where `n_points` is the
local volume and the 6 plaquette directions are:
`xy, xz, xt, yz, yt, zzt`. Each plaquette center is offset by 0.5 in
the two directions defining the plaquette.

### `get_extreme_plaq_xg_list`

```python
get_extreme_plaq_xg_list(
    plaq_xg_arr: np.ndarray,
    f_plaq: FieldRealD,
    threshold: float,
) -> list[tuple[float, list[float]]]
```

Return a globally sorted list of `(plaq, xg)` tuples for all plaquettes
with `plaq < 1 - threshold`. These are plaquette sites potentially
associated with instanton structures.

### `get_group_extreme_plaq_xg_list`

```python
get_group_extreme_plaq_xg_list(
    extreme_plaq_xg_list: list,
    total_site: Coordinate,
    dis_sqr_limit: float,
) -> list[list[tuple[float, list[float], float]]]
```

Group extreme plaquette points by proximity. Points within squared
distance `dis_sqr_limit` (using smearing-aware periodic distance) are
placed in the same group. Each group is represented as a list of
`(plaq, xg, dis_sq)` tuples, where `dis_sq` is the squared distance
from the group's first point.

---

## Instanton List Processing

### `process_inst_list`

```python
process_inst_list(inst_list: list) -> list[dict]
```

Process the raw instanton tracking data from `InstantonMap.inst_list`
into a structured list of per-instanton summary dictionaries.

Each entry in the returned `p_inst_list` contains:

| Key | Description |
|---|---|
| `inst_idx` | Original instanton index |
| `inst_raw` | Raw tracking data |
| `o_xg_d` | Coordinate on initial lattice (float) |
| `o_total_site` | Total site size on initial lattice |
| `xg_d` | Coordinate at detection (float) |
| `total_site` | Total site size at detection |
| `current_spacing` | Lattice spacing ratio |
| `plaq` | Minimum plaquette value |
| `flow_time` | Flow time at minimum plaq |
| `num_plaq` | Number of plaquettes in the group |
| `dis_sqr_max` | Max squared distance within group |
| `flow_time_start` / `flow_time_end` | Flow time range of tracking |
| `delta_s_topo` | Change in topological charge over flow |
| `estimate_topo_charge` | Estimated topological charge |
| `closest_inst_info` | Info about the closest instanton |
| `closest_sim_inst_info` | Info about the closest instanton with overlapping flow time |

### `get_tot_topo_count`

```python
get_tot_topo_count(p_inst_list: list[dict], flow_time: float) -> int
```

Get the total topological charge by summing the rounded estimated
topological charges of instantons detected after `flow_time`.

### `get_tot_inst_count`

```python
get_tot_inst_count(p_inst_list: list[dict], flow_time: float) -> int
```

Get the total number of instantons plus anti-instantons detected after
`flow_time`.

---

## Display Utilities

### `displayln_info_topo_info`

```python
displayln_info_topo_info(topo_info: dict) -> None
```

Print a formatted summary of a topological info record (flow time,
plaquette, topological charge, number of instantons, etc.).

### `displayln_info_p_inst`

```python
displayln_info_p_inst(p_inst: dict) -> None
```

Print a detailed summary of a single processed instanton record.

### `displayln_info_p_inst_list`

```python
displayln_info_p_inst_list(p_inst_list: list[dict]) -> None
```

Print summaries for all processed instantons and the cumulative
topological charge.

---

## InstantonMap Class

### `InstantonMap.__init__`

```python
InstantonMap(
    *,
    total_site: Coordinate = None,
    topo_sphere_sum_radius_list: list[float] = None,
    dis_sqr_limit: float = None,
    threshold: float = None,
)
```

Initialize an instanton tracker.

| Parameter | Description |
|---|---|
| `total_site` | Lattice dimensions (e.g., `q.Coordinate([4, 4, 4, 8])`) |
| `topo_sphere_sum_radius_list` | Radii for sphere-sum topological charge |
| `dis_sqr_limit` | Squared distance limit to identify same instanton |
| `threshold` | Plaquette threshold: sites with `plaq < 1 - threshold` are instanton candidates |

### `InstantonMap.shift`

```python
InstantonMap.shift(shift: Coordinate) -> None
```

Record a shift of the gauge field. Updates `origin_coordinate` so that
instanton coordinates remain in the initial lattice frame.

### `InstantonMap.convert_xg`

```python
InstantonMap.convert_xg(xg: Coordinate | CoordinateD) -> Coordinate | CoordinateD
```

Convert a coordinate in the current lattice frame to the initial lattice
frame.

### `InstantonMap.half_lattice`

```python
InstantonMap.half_lattice() -> None
```

Double the lattice spacing (halve the lattice size). Used when the gauge
field is coarsened during multi-scale instanton detection.

### `InstantonMap.acc_time`

```python
InstantonMap.acc_time(step_size: float) -> None
```

Accumulate flow time and increment the step counter after a flow step.

### `InstantonMap.acc_topo_info`

```python
InstantonMap.acc_topo_info(
    f_plaq: FieldRealD,
    f_topo: FieldRealD | None = None,
) -> None
```

Detect and track instantons from the current plaquette field. Groups
extreme plaquette sites, matches them to previously tracked instantons,
and updates the active instanton list. If `f_topo` is provided, computes
sphere-sum topological charge at each instanton location.

---

## Top-Level Computation

### `compute_inst_map`

```python
compute_inst_map(
    gf: GaugeField,
    *,
    dis_sqr_limit: float = None,
    threshold: float = None,
    topo_sphere_sum_radius_list: list[float] = None,
    plaq_min_threshold: float = None,
    max_spacing: int = None,
    flow_step_size_list: list[float] = None,
    flow_num_step_list: list[int] = None,
) -> dict
```

Run the full instanton map computation on a gauge field. Returns a
dictionary with keys `inst_list`, `info_list`, `flow_time_list`,
`tot_topo_list`.

The algorithm proceeds in three phases:
1. **Wilson flow** — smooth the field with standard Wilson flow
2. **Freeze flow** — freeze topological structure by suppressing plaquette deviations
3. **Shrink flow** — shrink instantons to detect and count them, optionally halving the lattice when the field is smooth enough

| Parameter | Default | Description |
|---|---|---|
| `dis_sqr_limit` | `5` | Squared distance for instanton grouping |
| `threshold` | `0.1` | Plaquette deviation threshold |
| `topo_sphere_sum_radius_list` | `[2, 3, 4, 5, 6]` | Radii for sphere-sum topological charge |
| `plaq_min_threshold` | `0.97^(1/16)` | Min plaq threshold for lattice halving |
| `max_spacing` | auto | Maximum lattice spacing (power of 2) |
| `flow_step_size_list` | `[0.05, 0.1, 0.1]` | Step sizes for Wilson, Freeze, Shrink |
| `flow_num_step_list` | `[80, 320, 1000000]` | Max steps for Wilson, Freeze, Shrink |

### `displayln_info_inst_map_obj`

```python
displayln_info_inst_map_obj(inst_map_obj: dict) -> None
```

Display a comprehensive summary of the instanton map results, including
per-instanton details and topological charge verification.

### `smear_measure_topo`

```python
smear_measure_topo(
    gf: GaugeField,
    smear_info_list: list[list] = None,
    *,
    energy_derivative_info: list = None,
    info_path: str = None,
    density_field_path: str = None,
    is_show_topo_terms: bool = False,
) -> tuple[list[dict], list[dict]]
```

Apply a sequence of smearing steps to `gf` in place, measuring
topological observables after each step. Returns `(topo_list, energy_list)`.

Each entry in `smear_info_list` is `[step_size, n_step, flow_type, wilson_flow_integrator_type]`.

The `topo_list` entries contain: `flow_time`, `plaq`, `plaq_min`,
`plaq_max`, `energy_density`, `energy_deriv`, `topo`, `topo_clf`,
`abs_topo`, and their per-timeslice variants.

The `energy_list` entries contain: `flow_time`, `plaq`, `plaq_min`,
`plaq_max`, `energy_density`, and per-timeslice values.

---

## Examples

### Compute Instanton Map

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

obj = q.compute_inst_map(gf)

inst_list = obj["inst_list"]
info_list = obj["info_list"]
flow_time_list = obj["flow_time_list"]
tot_topo_list = obj["tot_topo_list"]

print(f"Number of detected instantons: {len(inst_list)}")
print(f"Topological charge history: {tot_topo_list}")

q.end_with_mpi()
```

### Process and Display Instanton List

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

obj = q.compute_inst_map(gf)
q.displayln_info_inst_map_obj(obj)

q.end_with_mpi()
```

### Smear and Measure Topology

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

smear_info_list = [
    [0.05, 20, "Wilson", "euler"],
    [0.05, 20, "Wilson", "euler"],
    [0.01, 50, "DBW2", "euler"],
]

topo_list, energy_list = q.smear_measure_topo(gf, smear_info_list)

for entry in topo_list:
    print(f"t={entry['flow_time']:.4f}  topo={entry['topo']:.4f}  plaq={entry['plaq']:.6f}")

q.end_with_mpi()
```

### Plaq-Dependent Flow with Different Flow Types

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

step_size = 0.05

# Wilson flow
gf_w = gf.copy()
q.gf_flow_topo(gf_w, step_size, "Wilson")
print(f"After Wilson flow: plaq={gf_w.plaq():.6f}")

# Freeze flow
gf_f = gf.copy()
q.gf_flow_topo(gf_f, step_size, "Freeze")
print(f"After Freeze flow: plaq={gf_f.plaq():.6f}")

# Shrink flow
gf_s = gf.copy()
q.gf_flow_topo(gf_s, step_size, "Shrink")
print(f"After Shrink flow: plaq={gf_s.plaq():.6f}")

q.end_with_mpi()
```
