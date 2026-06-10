# `qlat.smear` — Gauge-Field and Propagator Smearing

Source: `qlat/qlat/smear.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Gauge-Field Smearing](#gauge-field-smearing)
   - [`gf_ape_smear`](#gf_ape_smear)
   - [`gf_spatial_ape_smear`](#gf_spatial_ape_smear)
   - [`gf_hyp_smear`](#gf_hyp_smear)
3. [Propagator and Fermion-Field Smearing](#propagator-and-fermion-field-smearing)
   - [`prop_smear`](#prop_smear)
   - [`prop_spatial_smear_no_comm`](#prop_spatial_smear_no_comm)
4. [Examples](#examples)

---

## Overview

The `qlat.smear` module provides smearing operations used in lattice QCD to
reduce UV fluctuations in gauge fields and to construct extended operators
for hadron spectroscopy. It wraps the C++ functions from `qcd-smear.h`.

Two families of operations are provided:

- **Gauge-field smearing** — APE and HYP smearing replace gauge links with
  weighted averages of surrounding paths, producing smoother gauge
  configurations suitable for computing topological charge, or for use as
  input to propagator smearing.
- **Propagator / fermion-field smearing** — Jacobi-style smearing applied to
  quark propagators and fermion fields, used to construct extended hadron
  operators with improved overlap onto ground-state hadrons.

All functions return new objects and leave the inputs unchanged.

---

## Gauge-Field Smearing

### `gf_ape_smear`

```python
gf_ape_smear(gf: GaugeField, alpha: float, steps: int = 1) -> GaugeField
```

Perform APE (Average Plaquette Evaluation) smearing on a gauge field.

Each step replaces each link U_mu(x) with a weighted combination of the
original link and the staples surrounding it:

```
U'_mu(x) = (1 - alpha) * U_mu(x) + (alpha / 6) * sum_of_staples
```

The result is re-projected to SU(3). After `steps` iterations the smeared
field is returned.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field |
| `alpha` | `float` | — | Smearing weight (0 = no smearing, 1 = full staple) |
| `steps` | `int` | `1` | Number of smearing iterations |
| **Returns** | `GaugeField` | | Smeared gauge field (new object) |

---

### `gf_spatial_ape_smear`

```python
gf_spatial_ape_smear(gf: GaugeField, alpha: float, steps: int = 1) -> GaugeField
```

Perform spatial-only APE smearing. Same as `gf_ape_smear` but only spatial
staples (mu, nu in {0, 1, 2}) are used; temporal links are unchanged.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field |
| `alpha` | `float` | — | Smearing weight |
| `steps` | `int` | `1` | Number of smearing iterations |
| **Returns** | `GaugeField` | | Smeared gauge field (new object) |

Typical values: `alpha = 0.5`, `steps = 30`.

---

### `gf_hyp_smear`

```python
gf_hyp_smear(gf: GaugeField, alpha1: float, alpha2: float, alpha3: float) -> GaugeField
```

Perform HYP (Hypersmearing) smearing. HYP smearing uses three nested levels
of staple construction with independent coefficients, providing stronger UV
suppression with less link distortion than APE smearing.

Reference: [Hasenfratz & Knechtli, Phys. Rev. D 64, 034504 (2001)](https://doi.org/10.1103/PhysRevD.64.034504), Eq. (4).

| Parameter | Type | Description |
|---|---|---|
| `gf` | `GaugeField` | Input gauge field |
| `alpha1` | `float` | Outer-level smearing coefficient |
| `alpha2` | `float` | Middle-level smearing coefficient |
| `alpha3` | `float` | Inner-level smearing coefficient |
| **Returns** | `GaugeField` | | Smeared gauge field (new object) |

Values from the paper: `alpha1 = 0.75`, `alpha2 = 0.6`, `alpha3 = 0.3`.

---

## Propagator and Fermion-Field Smearing

### `prop_smear`

```python
prop_smear(
    prop: Prop,
    gf1: GaugeField,
    coef: float,
    step: int,
    mom: CoordinateD = None,
    smear_in_time_dir: bool = False,
    mode_smear: int = 1,
) -> Prop
```

Apply Jacobi smearing to a propagator. The smearing operation at each
iteration is:

```
psi'(x) = coef * psi(x) + (1 - coef) / (2*3) * sum_{mu,nu} [U_mu(x) psi(x+mu) + U_mu^dag(x-mu) psi(x-mu)]
```

A momentum-dependent phase `exp(i * mom * x)` can be applied for
momentum-projected smearing.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `prop` | `Prop` | — | Input propagator |
| `gf1` | `GaugeField` | — | Gauge field (should be left-expanded, see below) |
| `coef` | `float` | — | Smearing coefficient (0 = full smearing, 1 = no smearing) |
| `step` | `int` | — | Number of smearing iterations |
| `mom` | `CoordinateD` or `None` | `None` | Momentum in lattice units 1/a (not 2*pi/L/a). `None` = zero momentum |
| `smear_in_time_dir` | `bool` | `False` | If `True`, include temporal direction in smearing |
| `mode_smear` | `int` | `1` | `0` = legacy `prop_smear`; `>= 1` = `prop_smear_qlat_convension` |
| **Returns** | `Prop` | | Smeared propagator (new object) |

**Typical workflow:**

```python
gf_ape = q.gf_spatial_ape_smear(gf, 0.5, 30)
gf1 = q.mk_left_expanded_field(gf_ape)
prop_smeared = q.prop_smear(prop, gf1, 0.9375, 10)
```

**Typical parameters by ensemble:**

| Ensemble | `coef` | `step` | Notes |
|---|---|---|---|
| 24D | `0.9375` | `10` | |
| 48I | `0.9375` | `29` | |

Set `mom = 0.5 * hadron_momentum` for momentum-projected smearing.

---

### `prop_spatial_smear_no_comm`

```python
prop_spatial_smear_no_comm(
    ff_list: list[FermionField4d],
    gf: GaugeField,
    coef: float,
    step: int,
    mom: CoordinateD = None,
) -> None
```

Apply spatial Jacobi smearing to a list of `FermionField4d` objects **in
place** without any MPI communication. Both `gf` and each field in
`ff_list` must contain entire time slices (no halo exchange is performed).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ff_list` | `list[FermionField4d]` | — | List of fermion fields to smear in place |
| `gf` | `GaugeField` | — | Gauge field covering the same geometry |
| `coef` | `float` | — | Smearing coefficient |
| `step` | `int` | — | Number of smearing iterations |
| `mom` | `CoordinateD` or `None` | `None` | Momentum in lattice units. `None` = zero momentum |

This function modifies `ff_list` in place and returns `None`. It is useful
when all time slices are local and the overhead of halo exchange is
unnecessary.

---

## Examples

### APE and HYP Smearing

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Load or generate a gauge field
gf = q.GaugeField(geo)
gf.set_unit()  # placeholder

# APE smearing: alpha=0.5, 30 iterations
gf_ape = q.gf_spatial_ape_smear(gf, 0.5, 30)
print(f"APE-smeared field: {gf_ape.crc32()}")

# HYP smearing with paper parameters
gf_hyp = q.gf_hyp_smear(gf, 0.75, 0.6, 0.3)
print(f"HYP-smeared field: {gf_hyp.crc32()}")

q.end_with_mpi()
```

### Propagator Smearing for Hadron Spectroscopy

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Prepare gauge field
gf = q.GaugeField(geo)
gf.set_unit()

# Smear gauge field spatially
gf_ape = q.gf_spatial_ape_smear(gf, 0.5, 30)
gf1 = q.mk_left_expanded_field(gf_ape)

# Smear a propagator
prop = q.Prop(geo)
prop.set_zero()

# 24D-like parameters, zero momentum
prop_smeared = q.prop_smear(prop, gf1, 0.9375, 10)

# With momentum projection (half the hadron momentum)
mom = q.CoordinateD([0.5, 0.0, 0.0, 0.0])
prop_mom_smeared = q.prop_smear(prop, gf1, 0.9375, 10, mom=mom)

q.end_with_mpi()
```

### In-Place Fermion Field Smearing

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

# prop_spatial_smear_no_comm requires size_node == (1, 1, 1, t_size),
# i.e., each MPI node must hold entire time slices. With a single node
# and size_node = (1, 1, 1, 1), t_size must be 1.
total_site = q.Coordinate([4, 4, 4, 1])
geo = q.Geometry(total_site)

gf = q.GaugeField(geo)
gf.set_unit()

# Create multiple fermion fields
ff_list = [q.FermionField4d(geo) for _ in range(3)]
for ff in ff_list:
    ff.set_zero()

# Smear all in place (no MPI communication)
q.prop_spatial_smear_no_comm(ff_list, gf, 0.9375, 10)

q.end_with_mpi()
```
