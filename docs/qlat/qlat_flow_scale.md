# `qlat.flow_scale` — Flow and Lattice Scale Determination

Source: `qlat/qlat/flow_scale.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Plaq Factor Construction](#plaq-factor-construction)
   - [`get_plaq_factor_for_gf_scale_flow`](#get_plaq_factor_for_gf_scale_flow)
3. [Scale Flow](#scale-flow)
   - [`gf_flow_scale`](#gf_flow_scale)
4. [Time-Slice Observables](#time-slice-observables)
   - [`gf_plaq_tslice`](#gf_plaq_tslice)
   - [`gf_energy_density_dir_tslice`](#gf_energy_density_dir_tslice)
5. [Flow Recording and Running](#flow-recording-and-running)
   - [`gf_flow_record`](#gf_flow_record)
   - [`run_flow_scale`](#run_flow_scale)
   - [`default_run_flow_scale_params`](#default_run_flow_scale_params)
6. [Examples](#examples)

---

## Overview

`flow_scale` generalizes the standard Wilson flow to support spatial-only
flow and defines additional flow observables for lattice scale
determination.

The standard Wilson flow evolves all plaquette directions uniformly. This
module allows restricting the flow to spatial plaquettes only (excluding
those involving the time direction), which is useful for scale setting
procedures that require anisotropic flow.

The module defines three flow observables:
$$
W_0(t) = t^2 \langle E(t) \rangle
$$
$$
W_1(t) = W(t) = t \frac{d}{dt} (t^2 \langle E(t) \rangle)
$$
$$
W_2(t) = -t^3 \frac{d}{dt} \langle E(t) \rangle
$$

Key references:
- [arXiv:1006.4518](https://arxiv.org/abs/1006.4518) — Wilson flow and energy density
- [arXiv:1203.4469](https://arxiv.org/abs/1203.4469) — Topological susceptibility from Wilson flow

---

## Plaq Factor Construction

### `get_plaq_factor_for_gf_scale_flow`

```python
get_plaq_factor_for_gf_scale_flow(
    total_site: tuple,
    is_spatial: bool,
    t_dir: int,
) -> FieldRealD
```

Construct a plaquette factor field for use with `gf_plaq_flow_force`.
Returns a `FieldRealD` with multiplicity 6 (one factor per plaquette
direction). The result is cached by `(total_site, is_spatial, t_dir)`.

When `is_spatial=False`, all factors are set to 1 (standard Wilson flow).
When `is_spatial=True`, only spatial plaquettes involving `t_dir` are
included (factors set to 1), and plaquettes involving the flow direction
are excluded (factors set to 0).

The 6 plaquette directions are ordered: `xy, xz, xt, yz, yt, zt`.

| `t_dir` | Spatial flow excludes | Included directions |
|---|---|---|
| 3 (time) | xt, yt, zt | xy, xz, yz |
| 0 (z) | xy, xz, xt | yz, yt, zt |
| 1 (y) | xy, yz, yt | xz, xt, zt |
| 2 (x) | xz, yz, zt | xy, xt, yt |

---

## Scale Flow

### `gf_flow_scale`

```python
gf_flow_scale(
    gf: GaugeField,
    step_size: float,
    *,
    is_spatial: bool = None,
    t_dir: int = None,
    integrator_type: str = None,
) -> None
```

Perform one step of scale flow on `gf` in place. Defaults to standard
(isotropic) Wilson flow with the 3rd-order Runge-Kutta integrator.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Gauge field (modified in place) |
| `step_size` | `float` | — | Flow step size |
| `is_spatial` | `bool` | `False` | If `True`, flow only spatial plaquettes |
| `t_dir` | `int` | `3` | Direction treated as "time" (0=x, 1=y, 2=z, 3=t) |
| `integrator_type` | `str` | `"runge-kutta"` | `"euler"` or `"runge-kutta"` |

The Runge-Kutta scheme follows [arXiv:1006.4518v3](http://arxiv.org/abs/1006.4518v3).

---

## Time-Slice Observables

### `gf_plaq_tslice`

```python
gf_plaq_tslice(gf: GaugeField, *, t_dir: int = None) -> np.ndarray
```

Compute the average plaquette per time-slice. Returns an array of shape
`(t_size, 6)` where `t_size = geo.total_site[t_dir]` and the 6 values
correspond to the plaquette directions `xy, xz, xt, yz, yt, zt`.

### `gf_energy_density_dir_tslice`

```python
gf_energy_density_dir_tslice(gf: GaugeField, *, t_dir: int = None) -> np.ndarray
```

Compute the average energy density per time-slice, decomposed by
direction. Returns an array of shape `(t_size, 6)`.

The relation to the plaquette is approximately:
```
gf_energy_density_dir_tslice(gf)[:]  ~  6 * (1 - gf_plaq_tslice(gf)[:])
```

---

## Flow Recording and Running

### `gf_flow_record`

```python
gf_flow_record(
    gf: GaugeField,
    step_size: float,
    num_step: int,
    *,
    flow_time: float = None,
    is_spatial: bool = None,
    t_dir: int = None,
    integrator_type: str = None,
) -> dict
```

Run `num_step` flow steps on `gf` in place, recording time-slice
observables after each step. Returns a dictionary:

```python
obj_record = dict(
    info_list=[
        dict(flow_time, plaq_tslice, energy_density_dir_tslice),
        ...
    ],
    params=dict(step_size, num_step, flow_time, is_spatial, t_dir, integrator_type),
)
```

The gauge field is unitarized after each step.

### `run_flow_scale`

```python
run_flow_scale(
    fn_out: str,
    *,
    get_gf: callable = None,
    fn_gf: str = None,
    params: dict = None,
) -> None
```

High-level driver that loads a gauge field, runs `gf_flow_record`, and
saves the result to a pickle file. Skips if the output file already
exists.

Either `get_gf` (a callable returning a `GaugeField`) or `fn_gf` (a
file path to load) must be provided, but not both.

| Parameter | Description |
|---|---|
| `fn_out` | Output pickle file path (must end with `.pickle`) |
| `get_gf` | Callable returning a `GaugeField` |
| `fn_gf` | Path to a gauge field file |
| `params` | Override parameters (merged with `default_run_flow_scale_params`) |

### `default_run_flow_scale_params`

```python
default_run_flow_scale_params = dict(
    step_size=0.05,
    num_step=400,
    is_spatial=False,
    t_dir=3,
    integrator_type="runge-kutta",
)
```

Default parameters for `run_flow_scale`.

---

## Examples

### Standard Wilson Flow with Scale Recording

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
num_step = 100

obj_record = q.gf_flow_record(gf, step_size, num_step)

for info in obj_record["info_list"]:
    t = info["flow_time"]
    plaq_ts = info["plaq_tslice"]
    print(f"t={t:.4f}  plaq_avg={plaq_ts.mean():.6f}")

q.end_with_mpi()
```

### Spatial-Only Flow

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

# Flow only spatial plaquettes (excluding time direction)
q.gf_flow_scale(gf, 0.05, is_spatial=True, t_dir=3, integrator_type="runge-kutta")

print(f"Plaq after spatial flow: {gf.plaq():.6f}")

q.end_with_mpi()
```

### Time-Slice Observables

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

plaq_ts = q.gf_plaq_tslice(gf)
ed_ts = q.gf_energy_density_dir_tslice(gf)

print(f"Plaq timeslice shape: {plaq_ts.shape}")
print(f"Energy density timeslice shape: {ed_ts.shape}")

for t in range(plaq_ts.shape[0]):
    print(f"t={t}  plaq_sum={plaq_ts[t].sum():.6f}  ed_sum={ed_ts[t].sum():.6f}")

q.end_with_mpi()
```

### Run Flow Scale with Custom Parameters

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

params = dict(
    step_size=0.02,
    num_step=200,
    is_spatial=False,
    t_dir=3,
    integrator_type="runge-kutta",
)

def get_gf():
    return gf.copy()

q.run_flow_scale("/tmp/flow-scale.pickle", get_gf=get_gf, params=params)

q.end_with_mpi()
```
