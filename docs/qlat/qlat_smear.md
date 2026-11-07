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
- The iteration count is `num_step` for the projection-based gauge-field
  smearings and for the propagator smearings. It is the last positional
  parameter and defaults to `1`; `num_step == 0` returns an unchanged copy
  (negative values behave as `0`). `gf_ape_smear_no_proj` is the exception:
  it takes a target accumulated alpha instead (see
  [`gf_ape_smear_no_proj`](#gf_ape_smear_no_proj)).
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
    gf: GaugeField, alpha: float, accumulated_alpha: float
) -> GaugeField
```

Perform an APE-like smearing that omits the SU(3) projection step. Unlike the
other gauge-field smearings this one is not driven by an iteration count:
`accumulated_alpha` is a **target** total effective smearing, and the function
performs as many steps as needed to accumulate that total (or until it can
make no further progress, see below).

Each step adds the staple term to the current link, without damping the link
term and without projecting the result back to SU(3), using a step-dependent
coefficient `a_k`:

```
U_{k+1, mu}(x) = U_{k, mu}(x) + (a_k / 6) * staple_{k, mu}(x)
```

This is not the convex interpolation used by `gf_ape_smear`: the linear term is
not scaled down, and for a near-unitary field the staple satisfies
`staple_{k, mu}(x) ~ 6 U_{k, mu}(x)`, so one step is approximately
`U_{k+1} ~ (1 + a_k) U_k` and the amplitude of the field **grows** with every
step.

The staple `staple_{k, mu}(x)` is built from the current field `U_k`, and

```
M_k      = max_{x, mu} Tr( U_{k, mu}(x)^dag U_{k, mu}(x) )
factor_k = 3 / M_k = 1 / max( Tr / 3 )
a_k      = min( alpha * factor_k, accumulated_alpha - sum_{j < k} a_j )
```

`M_k` is a **global** (MPI-reduced) maximum over all links, so every node
obtains the same `factor_k`. Nothing normalises the field itself: there is no
per-step or final multiplication by a common scalar, and the amplitude enters
the iteration only through `factor_k`, which sets the step size relative to
the current field magnitude. A unitary field has `M_0 = 3` and `factor_0 = 1`,
so the first step uses `a_0 = min(alpha, accumulated_alpha)`, an ordinary
unprojected step of weight `alpha` whenever the target is at least `alpha`.
For a near-unitary input that update grows the amplitude, so `M_k` rises above
`3`, `factor_k = 3 / M_k` falls below `1`, and every later step has
`a_k <= alpha`: the step size is only ever reduced, never enlarged, and it
cannot reach or exceed `1`.

`sum_{j < k} a_j` is the coefficient accumulated by the steps already applied,
and the last step is truncated to the remaining
`accumulated_alpha - sum_{j < k} a_j`; the loop therefore never overshoots the
target, and it stops as soon as the sum has reached it (up to floating-point
rounding, handled by the progress guard below). `accumulated_alpha` thus
measures the total weight applied, not the number of updates:
`accumulated_alpha = alpha * n` asks for the same total as `n` nominal steps,
but it is not reached in `n` updates. Because the field grows, `a_k` falls as
the amplitude rises and the accumulated total increases by much less than
`alpha` per step, so a given target can need far more steps than
`accumulated_alpha / alpha`: with `alpha = 0.5`, a target of `2.5` is
approached with on the order of a hundred small steps rather than five. Since
nothing renormalises the field, its growth can also make the sum of the `a_k`
converge, capping the accumulated weight below the target; the loop then stops
through the progress guard below, with the amplitude already grown by many
orders of magnitude.

The input `gf` is not modified; the result is returned as a new object. If
`alpha <= 0` or `accumulated_alpha <= 0` the loop exits before the first
update, so the returned copy is bit-exactly the input (in particular
`alpha = 0` and `accumulated_alpha = 0` are identities).

The loop also stops early when the field degenerates (`M_k <= 0` or not
finite, so `factor_k` cannot be formed) or when floating-point arithmetic can
no longer make progress towards the target -- the accumulated sum stops
increasing, which is what happens when the field growth caps the accumulated
weight below the target. At verbose level 0 and above (`q_verbose >= 0`), the
progress of the iteration -- the step index, `a_k`, the accumulated total and
`M_k` -- is printed every fourth step and on the final step. When the loop ends
below the target, the shortfall is always printed on screen, because it is a
level `-1` message: it appears even at the library default verbosity
(`q_verbose` unset, i.e. `-1`). The message reads
`gf_ape_smear_no_proj: reached only <acc> of accumulated_alpha = <target>`.
Neither message carries a `CHECK:`, `INFO:` or `WARNING` marker, so neither
reaches the harvested `<name>.log` / `.check.txt` reference files.

> **Warning:** the update is additive, so the field amplitude is not
> conserved: for near-unitary links `staple_k ~ 6 U_k`, and one step is
> roughly `U_{k+1} ~ (1 + a_k) U_k`. Nothing in this function rescales or
> renormalises the field, so the links grow away from SU(3) and the returned
> field is far from unitary. Because the field grows, `M_k > 3` for a
> near-unitary input, `factor_k = 3 / M_k < 1` and `a_k <= alpha`: unlike the
> interpolation in `gf_ape_smear`, the step size is then bounded by `alpha`,
> and the link term always keeps its full weight of `1`. The price is paid in
> the step count and in the reachable range: `a_k` decreases as the amplitude
> grows, the applied weight converges, and a target above the resulting cap
> is never reached -- a target such as `2.5` at `alpha = 0.5` can cost on the
> order of a hundred tiny steps and still fall short, leaving a finite but
> enormously large amplitude.
> `accumulated_alpha` measures the total weight applied along this
> unprojected, non-unitary trajectory, not a number of steps and not a
> property of a projected field; project the result back to SU(3) afterwards if
> that is required.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf` | `GaugeField` | — | Input gauge field (not modified) |
| `alpha` | `float` | — | Nominal smearing weight per step, scaled by `factor_k = 3 / M_k`; since the field grows, `M_k > 3` and `a_k <= alpha` |
| `accumulated_alpha` | `float` | — | **Target** total effective smearing (not a step count); the applied `a_k` are summed until it is reached (last step truncated), or until the progress guard stops the loop |
| **Returns** | `GaugeField` | | Smeared gauge field (new object) |

In C++ the same operation writes its result into `gf` (starting from
`gf = gf0`, with `gf == gf0` aliasing allowed) and returns nothing; the target
is the fourth argument `accumulated_alpha`.

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

Reference: [Hasenfratz & Knechtli, Phys. Rev. D 64, 034504 (2001)](https://doi.org/10.1103/PhysRevD.64.034504), equation (4).

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

# APE-like smearing without the SU(3) projection, driven by a target total
# smearing instead of a step count. The step is additive, U' = U + (a/6) *
# staple, so the field amplitude grows and the weight per step decreases;
# alpha = 0.5 with accumulated_alpha = 1.5 targets the same total weight as
# three nominal steps of 0.5, approached through many small steps (the growth
# can cap the accumulated weight below the target).
gf_no_proj = q.gf_ape_smear_no_proj(gf, 0.5, 1.5)
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
