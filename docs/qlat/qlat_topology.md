# `qlat.topology` — Topology Charge and Plaquette Action Density

Source: `qlat/qlat/topology.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Plaquette Action Density](#plaquette-action-density)
   - [`gf_plaq_action_density_field`](#gf_plaq_action_density_field)
   - [`gf_plaq_action_density`](#gf_plaq_action_density)
   - [`gf_spatial_plaq_action_density_field`](#gf_spatial_plaq_action_density_field)
   - [`gf_spatial_plaq_action_density`](#gf_spatial_plaq_action_density)
3. [Topological Charge](#topological-charge)
   - [`gf_topology_field_clf`](#gf_topology_field_clf)
   - [`gf_topology_clf`](#gf_topology_clf)
   - [`gf_topology_field`](#gf_topology_field)
   - [`gf_topology`](#gf_topology)
   - [`gf_topology_terms_field`](#gf_topology_terms_field)
   - [`gf_topology_terms`](#gf_topology_terms)
4. [Examples](#examples)

---

## Overview

`topology` provides measurements of the plaquette action density and
topological charge on gauge field configurations. Two definitions of the
topological charge are available:

- **Clover-leaf** (`_clf` variants): uses the basic clover discretization of
  the field strength tensor.
- **5-loop improved** (default variants): uses the 5-loop improved definition
  from [arXiv:hep-lat/9701012](https://arxiv.org/abs/hep-lat/9701012),
  which reduces lattice artifacts.

The plaquette action density functions compute the sum over all six plaquettes
at each site of `(1 - 1/3 * Re Tr U_P)`, which is the standard Wilson action
density. The total action is `beta * total_volume * action_density`.

For a single instanton the action is `8 * pi^2 / g^2`, and `beta = 6 / g^2`.

---

## Plaquette Action Density

### `gf_plaq_action_density_field`

```python
gf_plaq_action_density_field(gf: GaugeField) -> FieldRealD
```

Return a `FieldRealD` with `multiplicity == 1` containing the plaquette
action density at each site:

```
paf(x) = sum_P (1 - 1/3 * Re Tr U_P)
```

where the sum runs over the 6 plaquettes originating at site `x`. The total
action is:

```
Action = beta * total_volume * paf_glb_sum / total_volume
```

### `gf_plaq_action_density`

```python
gf_plaq_action_density(gf: GaugeField) -> float
```

Return the spatially averaged plaquette action density as a scalar. Equivalent
to:

```python
gf_plaq_action_density_field(gf).glb_sum()[:].item() / total_volume
```

### `gf_spatial_plaq_action_density_field`

```python
gf_spatial_plaq_action_density_field(gf: GaugeField) -> FieldRealD
```

Return a `FieldRealD` with `multiplicity == 1` containing the spatial
plaquette action density at each site. Only spatial plaquettes (those not
involving the time direction) are included in the sum.

### `gf_spatial_plaq_action_density`

```python
gf_spatial_plaq_action_density(gf: GaugeField) -> float
```

Return the spatially averaged spatial-plaquette action density as a scalar.

---

## Topological Charge

### `gf_topology_field_clf`

```python
gf_topology_field_clf(gf: GaugeField) -> FieldRealD
```

Return a `FieldRealD` with `multiplicity == 1` containing the topological
charge density at each site using the basic clover-leaf discretization of
the field strength tensor. This is **not** the improved definition.

### `gf_topology_clf`

```python
gf_topology_clf(gf: GaugeField) -> float
```

Return the total topological charge as a scalar, using the basic clover-leaf
definition. Equivalent to:

```python
gf_topology_field_clf(gf).glb_sum()[:].item()
```

### `gf_topology_field`

```python
gf_topology_field(gf: GaugeField) -> FieldRealD
```

Return a `FieldRealD` with `multiplicity == 1` containing the topological
charge density at each site using the 5-loop improved definition from
[arXiv:hep-lat/9701012 Eq. (2-7)](https://arxiv.org/abs/hep-lat/9701012).
This definition reduces discretization errors compared to the basic
clover-leaf method.

### `gf_topology`

```python
gf_topology(gf: GaugeField) -> float
```

Return the total topological charge as a scalar, using the 5-loop improved
definition. Equivalent to:

```python
gf_topology_field(gf).glb_sum()[:].item()
```

### `gf_topology_terms_field`

```python
gf_topology_terms_field(gf: GaugeField) -> FieldRealD
```

Return a `FieldRealD` with `multiplicity == 5` containing the five individual
loop terms that sum to give the 5-loop improved topological charge density.
Useful for analyzing the contribution of each term. The sum over the last
axis equals `gf_topology_field`:

```python
gf_topology_terms_field(gf)[:].sum(-1) == gf_topology_field(gf)[:]
```

### `gf_topology_terms`

```python
gf_topology_terms(gf: GaugeField) -> np.ndarray
```

Return a NumPy array of shape `(5,)` and dtype `float64` containing the
globally summed values of each of the five loop terms. The sum of the array
equals `gf_topology(gf)`.

---

## Examples

### Plaquette Action Density

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

pa = q.gf_plaq_action_density(gf)
print(f"Plaquette action density: {pa}")

pa_field = q.gf_plaq_action_density_field(gf)
print(f"Per-site shape: {pa_field.geo.multiplicity}")

q.end_with_mpi()
```

### Topological Charge (5-loop improved)

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

top = q.gf_topology(gf)
print(f"Topological charge (5-loop): {top}")

q.end_with_mpi()
```

### Comparing Clover-Leaf and 5-loop Topology

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

top_clf = q.gf_topology_clf(gf)
top_5loop = q.gf_topology(gf)

print(f"Clover-leaf topology: {top_clf}")
print(f"5-loop improved:      {top_5loop}")

q.end_with_mpi()
```

### Decomposing the 5-loop Terms

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

terms = q.gf_topology_terms(gf)
print(f"5-loop terms: {terms}")
print(f"Sum of terms: {terms.sum()}")
print(f"Direct 5-loop: {q.gf_topology(gf)}")

q.end_with_mpi()
```

### Spatial Plaquette Action Density

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

pa_total = q.gf_plaq_action_density(gf)
pa_spatial = q.gf_spatial_plaq_action_density(gf)

print(f"Total action density:   {pa_total}")
print(f"Spatial action density: {pa_spatial}")

q.end_with_mpi()
```
