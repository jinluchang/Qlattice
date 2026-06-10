# `qlat.hlbl_contract` — Hadronic Light-by-Light Contractions

Source: `qlat/qlat/hlbl_contract.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Muon-Line Field](#muon-line-field)
3. [Pair Label Generators](#pair-label-generators)
4. [Local Current from Propagators](#local-current-from-propagators)
5. [Point-Selected Probability Fields](#point-selected-probability-fields)
6. [`CurrentMoments` Class](#currentmoments-class)
   - [Constructor](#constructor)
   - [Copy and Assignment](#copy-and-assignment)
   - [Compute Moments](#compute-moments)
   - [Global Sum](#global-sum)
7. [Four-Pair Contraction](#four-pair-contraction)
8. [Two-Plus-Two Pair Contraction](#two-plus-two-pair-contraction)
9. [Examples](#examples)

---

## Overview

`hlbl_contract` implements the hadronic light-by-light (HLbL) scattering
contraction routines needed for the muon anomalous magnetic moment (g-2)
calculation on the lattice. The module provides:

- **Muon-line fields** (`mk_m_z_field_tag`) — the electromagnetic vertex
  kernel at a point z, used to connect the muon line to the hadronic
  vacuum polarization.
- **Local currents** (`mk_local_current_from_props`) — the electromagnetic
  current J_mu(x) constructed from quark propagator pairs.
- **`CurrentMoments`** — moments of the current distribution used in the
  multipole expansion of the HLbL tensor.
- **Four-pair contractions** (`contract_four_pair_no_glb_sum`) — the
  dominant HLbL diagram where two photon propagators connect to a single
  quark loop.
- **Two-plus-two pair contractions** (`contract_two_plus_two_pair_no_glb_sum`)
  — the subdominant HLbL diagram with two separate quark loops connected
  by two photon propagators.

All contraction functions return results as NumPy arrays with shape
`(n_labels, s_limit, l_limit)` where `s` and `l` are angular momentum
quantum numbers in the multipole expansion. Global sums over MPI ranks
are not performed inside these functions; the caller must sum the returned
arrays.

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])
# ... (see examples below)
q.end_with_mpi()
```

---

## Muon-Line Field

### `mk_m_z_field_tag(psel_d, xg_x, xg_y, a, tag) -> SelectedPointsRealD`

Compute the muon-line field at vertex z for sampled source point x and
sink point y.

| Parameter | Type | Description |
|---|---|---|
| `psel_d` | `PointsSelection` | Point selection for the d (detector) sites |
| `xg_x` | `Coordinate` | Global coordinate of source point x |
| `xg_y` | `Coordinate` | Global coordinate of sink point y |
| `a` | `float` | Lattice spacing (in muon-mass units; use `muon_mass * a` in lattice units) |
| `tag` | `int` | `0` = subtraction (preferred), `1` = no subtraction |

Returns a `SelectedPointsRealD` containing the muon-line field values at
the selected points.

---

## Pair Label Generators

### `contract_four_pair_labels(tags: list[str]) -> list`

Generate the contraction labels for the four-pair HLbL diagram. Each tag
in `tags` specifies a reference point configuration:
- `"ref-far"` — reference point far from the interaction region
- `"ref-center"` — reference point at the center
- `"ref-close"` — reference point close to the interaction region

### `contract_two_plus_two_pair_labels() -> list`

Generate the contraction labels for the two-plus-two pair HLbL diagram.

---

## Local Current from Propagators

### `mk_local_current_from_props(sprop1, sprop2) -> SelectedPointsWilsonMatrix`

Construct the local electromagnetic current from a pair of point-selected
propagators:

```
J(x) = gamma5 * sprop2^+ * gamma5 * gamma_mu * sprop1
```

| Parameter | Type | Description |
|---|---|---|
| `sprop1` | `PselProp` | First propagator (forward) |
| `sprop2` | `PselProp` | Second propagator (conjugated with gamma5) |

Returns a `SelectedPointsWilsonMatrix` containing the current at each
selected point.

---

## Point-Selected Probability Fields

### `mk_psel_d_prob_xy(psel_prob, psel_d_prob, idx_xg_x, idx_xg_y) -> tuple`

Compute the probability-weighted field for a pair of source/sink points.

| Parameter | Type | Description |
|---|---|---|
| `psel_prob` | `SelectedPointsRealD` | Probability weights for the primary point selection |
| `psel_d_prob` | `SelectedPointsRealD` | Probability weights for the detector point selection |
| `idx_xg_x` | `int` | Index of source point x in `psel` |
| `idx_xg_y` | `int` | Index of sink point y in `psel` |

Returns `(prob_pair, psel_d_prob_xy)` where `prob_pair` is the combined
probability for the (x, y) pair and `psel_d_prob_xy` is the per-detector-point
probability field.

---

## `CurrentMoments` Class

`CurrentMoments` stores the multipole moments of the electromagnetic current
distribution, used in the HLbL tensor decomposition.

### Constructor

### `CurrentMoments()`

Create an uninitialized instance.

### `CurrentMoments(current: SelectedPointsWilsonMatrix, psel_d_prob_xy: SelectedPointsRealD)`

Create and compute moments from a local current field weighted by
detector-point probabilities.

---

### Copy and Assignment

### `copy(is_copying_data: bool = True) -> CurrentMoments`

Return a copy. If `is_copying_data` is `False`, return an uninitialized
instance.

### `__imatmul__(v1: CurrentMoments)`

In-place assignment: `self @= v1`.

---

### Compute Moments

### `set_from_current(current: SelectedPointsWilsonMatrix, psel_d_prob_xy: SelectedPointsRealD)`

Compute moments from the given current and probability fields.

---

### Global Sum

### `glb_sum() -> CurrentMoments`

Return a new `CurrentMoments` with all entries summed across MPI ranks.
Does not modify `self`.

---

## Four-Pair Contraction

### `contract_four_pair_no_glb_sum(coef, psel_prob, psel_d_prob, idx_xg_x, idx_xg_y, smf_d, sc_xy, sc_yx, cm_xy, cm_yx, inv_type, tags, r_sq_limit, muon_mass, z_v) -> numpy.ndarray`

Perform the four-pair HLbL contraction for a given (x, y) source/sink pair.
Returns a 3-D NumPy array of shape `(n_labels, s_limit, l_limit)`.

| Parameter | Type | Description |
|---|---|---|
| `coef` | `complex` | Overall coefficient (default 1.0) |
| `psel_prob` | `SelectedPointsRealD` | Probability weights for primary point selection |
| `psel_d_prob` | `SelectedPointsRealD` | Probability weights for detector point selection |
| `idx_xg_x` | `int` | Index of source point x in `psel` |
| `idx_xg_y` | `int` | Index of sink point y in `psel` |
| `smf_d` | `SelectedPointsRealD` | Muon-line field at vertex z |
| `sc_xy` | `SelectedPointsWilsonMatrix` | Local current J(x) for x->y direction |
| `sc_yx` | `SelectedPointsWilsonMatrix` | Local current J(x) for y->x direction |
| `cm_xy` | `CurrentMoments` | Current moments for x->y direction |
| `cm_yx` | `CurrentMoments` | Current moments for y->x direction |
| `inv_type` | `int` | `0` = light quark, `1` = strange quark |
| `tags` | `list[str]` | Reference point tags (e.g., `["ref-far", "ref-center"]`) |
| `r_sq_limit` | `int` | Maximum r^2 cutoff |
| `muon_mass` | `float` | Muon mass |
| `z_v` | `float` | Vector coupling constant |

**Note:** The returned array has not been globally summed; the caller
must sum across MPI ranks.

---

## Two-Plus-Two Pair Contraction

### `contract_two_plus_two_pair_no_glb_sum(coef, psel_prob, psel_lps_prob, idx_xg_x, lps_hvp_x, edl_list_c, r_sq_limit, muon_mass, z_v) -> tuple`

Perform the two-plus-two pair HLbL contraction. This diagram involves two
quark loops connected by two photon propagators.

Returns `(n_points_selected, n_points_computed, lsl_arr)` where `lsl_arr`
has shape `(n_labels, s_limit, l_limit)`.

| Parameter | Type | Description |
|---|---|---|
| `coef` | `complex` | Overall coefficient |
| `psel_prob` | `SelectedPointsRealD` | Probability weights for primary point selection |
| `psel_lps_prob` | `SelectedPointsRealD` | Probability weights for the loop point selection |
| `idx_xg_x` | `int` | Index of source point x in `psel` |
| `lps_hvp_x` | `SelectedPointsComplexD` | HVP tensor at point x |
| `edl_list_c` | `SelectedPointsComplexD` | External double-line list |
| `r_sq_limit` | `int` | Maximum r^2 cutoff |
| `muon_mass` | `float` | Muon mass |
| `z_v` | `float` | Vector coupling constant |

**Note:** The returned array has not been globally summed; the caller
must sum across MPI ranks.

---

## Examples

### Compute Muon-Line Field

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

psel = q.PointsSelection(geo, [q.Coordinate([0, 0, 0, 0])])
psel_d = q.PointsSelection(geo, [q.Coordinate([1, 0, 0, 0]),
                                  q.Coordinate([0, 1, 0, 0])])

xg_x = q.Coordinate([0, 0, 0, 0])
xg_y = q.Coordinate([2, 0, 0, 0])
a = 0.1  # lattice spacing in muon-mass units
tag = 0  # subtraction (preferred)

smf_d = q.hlbl_contract.mk_m_z_field_tag(psel_d, xg_x, xg_y, a, tag)
print(f"muon-line field shape: {smf_d.n_points}")

q.end_with_mpi()
```

### Local Current from Propagators

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

psel = q.PointsSelection(geo, [q.Coordinate([0, 0, 0, 0])])

sprop1 = q.PselProp(psel)
sprop2 = q.PselProp(psel)
sprop1.set_zero()
sprop2.set_zero()

scf = q.hlbl_contract.mk_local_current_from_props(sprop1, sprop2)
print(f"local current computed at {scf.n_points} points")

q.end_with_mpi()
```

### Four-Pair Contraction Workflow

```python
import qlat as q
import numpy as np

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Create minimal selections
psel = q.PointsSelection(geo, [q.Coordinate([0, 0, 0, 0]),
                                q.Coordinate([1, 0, 0, 0])])
psel_d = q.PointsSelection(geo, [q.Coordinate([0, 0, 0, 0])])

# Probability weights (uniform for demonstration)
psel_prob = q.SelectedPointsRealD(psel)
psel_d_prob = q.SelectedPointsRealD(psel_d)
psel_prob.set_zero()
psel_d_prob.set_zero()

# Generate pair labels
tags = ["ref-far", "ref-center"]
labels = q.hlbl_contract.contract_four_pair_labels(tags)
print(f"number of labels: {len(labels)}")

# Generate two-plus-two labels
labels_22 = q.hlbl_contract.contract_two_plus_two_pair_labels()
print(f"number of 2+2 labels: {len(labels_22)}")

q.end_with_mpi()
```

### CurrentMoments

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

psel = q.PointsSelection(geo, [q.Coordinate([0, 0, 0, 0])])
psel_d = q.PointsSelection(geo, [q.Coordinate([0, 0, 0, 0])])

# Build a current and probability field
sprop1 = q.PselProp(psel)
sprop2 = q.PselProp(psel)
sprop1.set_zero()
sprop2.set_zero()

scf = q.hlbl_contract.mk_local_current_from_props(sprop1, sprop2)
psel_d_prob = q.SelectedPointsRealD(psel_d)
psel_d_prob.set_zero()

cm = q.hlbl_contract.CurrentMoments(scf, psel_d_prob)
cm_glb = cm.glb_sum()
print("CurrentMoments computed and globally summed")

q.end_with_mpi()
```
