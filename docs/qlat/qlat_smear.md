# `qlat.smear` — Gauge-Field and Propagator Smearing

Source: `qlat/qlat/smear.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Gauge-Field Smearing](#gauge-field-smearing)
   - [`gf_ape_smear`](#gf_ape_smear)
   - [`gf_ape_smear_spatial`](#gf_ape_smear_spatial)
   - [`gf_ape_smear_no_proj`](#gf_ape_smear_no_proj)
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

All gauge-field smearing functions return new objects and leave the inputs
unchanged. `prop_smear` also returns a new object, while
`prop_spatial_smear_no_comm` smears the given fermion fields in place.

### Naming conventions

- Names are object-first: `<object>_<method>_smear[_<variant>]`, e.g.
  `gf_ape_smear_spatial`, `prop_smear`, `field_smear_step`.
- The iteration count is always `num_step`. It is the last positional
  parameter and defaults to `1`; `num_step == 0` returns an unchanged copy
  (negative values behave as `0`).
- Every public smearing function returns a new object and does not modify its
  input. The single exception is the internal low-level kernel
  `prop_spatial_smear_no_comm`, which works in place and returns `None`.

---

## Gauge-Field Smearing

### `gf_ape_smear`

```python
gf_ape_smear(gf: GaugeField, alpha: float, num_step: int = 1) -> GaugeField
```

Perform APE (Average Plaquette Evaluation) smearing on a gauge field.

Each step replaces each link U_mu(x) with a weighted combination of the
original link and the staples surrounding it:

```
U'_mu(x) = (1 - alpha) * U_mu(x) + (alpha / 6) * sum_of_staples
```

The result is re-projected to SU(3). After `num_step` iterations the smeared
field is returned.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field |
| `alpha` | `float` | — | Smearing weight (0 = no smearing, 1 = full staple) |
| `num_step` | `int` | `1` | Number of smearing iterations |
| **Returns** | `GaugeField` | | Smeared gauge field (new object) |

---

### `gf_ape_smear_no_proj`

```python
gf_ape_smear_no_proj(
    gf: GaugeField, alpha: float, num_step: int = 1
) -> GaugeField
```

Perform an APE-like smearing that omits the SU(3) projection step. Each
iteration uses the same weighted link/staple combination as `gf_ape_smear`:

```
U'_mu(x) = (1 - alpha) * U_mu(x) + (alpha / 6) * sum_of_staples
```

but the result is **not** re-projected to SU(3), so the update alone does not
pin the link amplitude. To compensate, after every iteration -- including the
last one -- the whole field is multiplied by one real scalar

```
s_k = sqrt( 3 / <N^2>_raw )

<N^2>_raw = (1 / (V * DIMN)) * sum_{x, mu} Tr( U_mu(x)^dag U_mu(x) )
```

where `V = geo.total_volume()`, `DIMN` is the number of spacetime dimensions
and `<N^2>_raw` is a **global** (MPI-reduced) average over all links, so every
node applies the same `s_k`. Step 0 is the raw input (`s_0 = 1`); after every
rescaled step this global average is `<N^2> = 3`, the value for a unitary link.
The input `gf` is not modified; with `num_step == 0` the returned field is an
unchanged copy.

The reported statistic is `N^2 = Tr(U^dag U) = sum_ab |U_ab|^2`, never its
square root: `N^2` is the natural link magnitude (exactly 3 for a unitary
link), it is the quantity the unprojected update shrinks, and it avoids taking
a square root of a quantity that can pass through zero when links vanish. The
square root appears only in `s_k`, because the correction acts on the amplitude
while the constraint is on the squared norm.

The rescale is unconditional: there is no parameter or flag that disables it.
In particular `alpha = 0` disables only the staple term, so a non-unitary input
is still rescaled towards `<N^2> = 3` at every step (for a unitary input
`s = 1` up to rounding). If the measured mean is not strictly positive (e.g. an
all-zero or NaN field), that step's rescale is skipped instead of producing
`inf`/`NaN`; this is the only case in which the field is not rescaled.

> **Warning:** the rescale sits inside the iteration, so it changes the
> update; it is not a pure post-processing factor. The map is not scale
> invariant (the linear term is first order in the field, the staple term
> third order), so removing the overall shrinkage keeps the cubic term at
> `O(1)` relative to the linear term and the iteration becomes genuinely
> nonlinear. At `alpha = 0.5` the field degenerates by roughly step 8: the
> unitarity deviation saturates at ~2.7, the mean `|det U|` falls to ~5e-6 at
> step 8 and ~1e-24 by step 16, and at step 32 most sites have `|det U| = 0`.
> There is no fixed scale factor either: measured values start at
> `s_1 = 1.224` (the first step grows the field slightly), fall to
> `s_4 = 1.0014`, then go below 1 from step 5 on (`s_k < 1`) and collapse to
> `s_8 ~ 2e-3`, while the cumulative product `S_k = s_1 * ... * s_k` is not
> monotone (it grows to ~1.359 at step 4, then decays). Do not use more than
> about 7 steps on this map at `alpha = 0.5`; the informative range is 1-4
> steps, and beyond that the configuration is not physical. If only a common
> overall amplitude is required, apply one factor after the whole ladder
> instead of running a long unprojected APE chain.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field |
| `alpha` | `float` | — | Smearing weight (0 = no staple mixing, 1 = full staple); the per-step rescale applies regardless |
| `num_step` | `int` | `1` | Number of smearing iterations |
| **Returns** | `GaugeField` | | Smeared gauge field (new object) |

---

### `gf_ape_smear_spatial`

```python
gf_ape_smear_spatial(
    gf: GaugeField, alpha: float, num_step: int = 1
) -> GaugeField
```

Perform spatial-only APE smearing. Same as `gf_ape_smear` but only spatial
staples (mu, nu in {0, 1, 2}) are used; temporal links are unchanged.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field |
| `alpha` | `float` | — | Smearing weight |
| `num_step` | `int` | `1` | Number of smearing iterations |
| **Returns** | `GaugeField` | | Smeared gauge field (new object) |

Typical values: `alpha = 0.5`, `num_step = 30`.

---

### `gf_hyp_smear`

```python
gf_hyp_smear(
    gf: GaugeField,
    alpha1: float,
    alpha2: float,
    alpha3: float,
    num_step: int = 1,
) -> GaugeField
```

Perform HYP (Hypersmearing) smearing. HYP smearing uses three nested levels
of staple construction with independent coefficients, providing stronger UV
suppression with less link distortion than APE smearing. After `num_step`
iterations the smeared field is returned.

Reference: [Hasenfratz & Knechtli, Phys. Rev. D 64, 034504 (2001)](https://doi.org/10.1103/PhysRevD.64.034504), Eq. (4).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field |
| `alpha1` | `float` | — | Outer-level smearing coefficient |
| `alpha2` | `float` | — | Middle-level smearing coefficient |
| `alpha3` | `float` | — | Inner-level smearing coefficient |
| `num_step` | `int` | `1` | Number of smearing iterations |
| **Returns** | `GaugeField` | | Smeared gauge field (new object) |

Values from the paper: `alpha1 = 0.75`, `alpha2 = 0.6`, `alpha3 = 0.3`.

---

## Propagator and Fermion-Field Smearing

### `prop_smear`

```python
prop_smear(
    prop: Prop,
    gf: GaugeField,
    coef: float,
    num_step: int = 1,
    mom: CoordinateD = None,
    smear_in_time_dir: bool = False,
    mode: str = "qlat",
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
| `gf` | `GaugeField` | — | Gauge field (should be left-expanded, see below) |
| `coef` | `float` | — | Smearing coefficient (0 = full smearing, 1 = no smearing) |
| `num_step` | `int` | `1` | Number of smearing iterations |
| `mom` | `CoordinateD` or `None` | `None` | Momentum in lattice units 1/a (not 2*pi/L/a). `None` = zero momentum |
| `smear_in_time_dir` | `bool` | `False` | If `True`, include temporal direction in smearing |
| `mode` | `str` | `"qlat"` | Smearing convention: `"qlat"` (default, qlat convention) or `"standard"` (legacy convention) |
| **Returns** | `Prop` | | Smeared propagator (new object) |

**Typical workflow:**

```python
gf_ape = q.gf_ape_smear_spatial(gf, 0.5, num_step=30)
gf_ext = q.mk_left_expanded_field(gf_ape)
prop_smeared = q.prop_smear(prop, gf_ext, 0.9375, num_step=10)
```

**Typical parameters by ensemble:**

| Ensemble | `coef` | `num_step` | Notes |
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

This is an internal, low-level kernel used by
`qlat.smear_prop.prop_spatial_smear`. It is deliberately excluded from the
unified smear API conventions: it keeps its original `step` argument, its
in-place semantics and its `None` return value. It is exported from `qlat.c`
for advanced use, but new code should normally call `prop_spatial_smear`
or `prop_smear` instead.

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
gf_ape = q.gf_ape_smear_spatial(gf, 0.5, num_step=30)
print(f"APE-smeared field: {gf_ape.crc32()}")

# APE-like smearing without the SU(3) projection. The field is rescaled to
# <N^2> = 3 after every step; keep the ladder short (a few steps at most).
gf_no_proj = q.gf_ape_smear_no_proj(gf, 0.5, num_step=3)
print(f"Unprojected APE-like field: {gf_no_proj.crc32()}")

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
gf_ape = q.gf_ape_smear_spatial(gf, 0.5, num_step=30)
gf_ext = q.mk_left_expanded_field(gf_ape)

# Smear a propagator
prop = q.Prop(geo)
prop.set_zero()

# 24D-like parameters, zero momentum
prop_smeared = q.prop_smear(prop, gf_ext, 0.9375, num_step=10)

# With momentum projection (half the hadron momentum)
mom = q.CoordinateD([0.5, 0.0, 0.0, 0.0])
prop_mom_smeared = q.prop_smear(prop, gf_ext, 0.9375, num_step=10, mom=mom)

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
