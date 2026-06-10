# `qlat.geometry` — Lattice Geometry and MPI Decomposition

Source: `qlat/qlat/geometry.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`Geometry` Class](#geometry-class)
   - [Constructors](#constructors)
   - [Properties](#properties)
   - [Coordinate Mapping](#coordinate-mapping)
   - [Site Queries](#site-queries)
   - [Display and Pickle](#display-and-pickle)
3. [`geo_resize` Function](#geo_resize-function)
4. [`geo_eo` Function](#geo_eo-function)
5. [Examples](#examples)

---

## Overview

`geometry` defines the `Geometry` class, which describes how a 4-D lattice is
decomposed across MPI nodes. A `Geometry` stores:

- The **total lattice dimensions** (`total_site = size_node * node_site`).
- The **local node site dimensions** (`node_site`) — the portion of the
  lattice owned by this MPI process.
- The **MPI grid layout** (`size_node`, `coor_node`, `id_node`).
- **Halo expansion** (`expansion_left`, `expansion_right`) — extra sites
  beyond the local volume needed for stencil operations.
- **Even/odd partitioning** (`eo`) — 0 for full lattice, 1 for odd sites,
  2 for even sites.

Every `Field` in qlat is defined on a `Geometry`. The geometry determines
which sites each MPI process stores, how indices map to coordinates, and
how halo regions are managed.

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
print(geo)                    # Geometry([4, 4, 4, 8])
print(geo.local_volume)      # 512
print(geo.total_volume)      # 512

q.end_with_mpi()
```

---

## `Geometry` Class

### Constructors

#### `Geometry()`

Create an uninitialized geometry.

#### `Geometry(total_site: Coordinate)`

Initialize from the total lattice dimensions. The MPI grid (`size_node`) and
this process's position (`coor_node`) are determined from the global `geon`
(set by `begin_with_mpi`).

| Parameter | Type | Description |
|---|---|---|
| `total_site` | `Coordinate` | Total lattice dimensions, e.g. `Coordinate([4, 4, 4, 8])` |

```python
geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))
```

#### `Geometry(id_node, size_node, node_site)`

Initialize from an integer node index without calling MPI functions.

| Parameter | Type | Description |
|---|---|---|
| `id_node` | `int` | Linear node index |
| `size_node` | `Coordinate` | MPI grid dimensions |
| `node_site` | `Coordinate` | Local site dimensions per node |

The relationship is:
```python
coor_node = q.coordinate_from_index(id_node, size_node)
id_node = q.index_from_coordinate(coor_node, size_node)
total_site = size_node * node_site
```

#### `Geometry(coor_node, size_node, node_site)`

Initialize from a coordinate node index without calling MPI functions.

| Parameter | Type | Description |
|---|---|---|
| `coor_node` | `Coordinate` | This node's position in the MPI grid |
| `size_node` | `Coordinate` | MPI grid dimensions |
| `node_site` | `Coordinate` | Local site dimensions per node |

---

### Properties

#### `total_site -> Coordinate`

The total lattice dimensions: `size_node * node_site`.

#### `total_volume -> int`

Total number of sites in the full lattice: `product(total_site)`.

#### `spatial_volume -> int`

Total number of spatial sites: `total_volume / total_site[3]`.

#### `local_site -> Coordinate`

Local site dimensions on this node. Same as `node_site`.

#### `local_volume -> int`

Number of sites on this node. For even/odd geometries (`eo != 0`), this is
half the product of `node_site`.

#### `local_volume_expanded -> int`

Number of sites including halo regions: `product(node_site_expanded)`.
For even/odd geometries, this is half.

#### `eo -> int`

Even/odd flag: `0` (full lattice), `1` (odd sites only), `2` (even sites
only).

#### `expansion_left -> Coordinate`

Halo expansion on the negative (left) side of each direction.

#### `expansion_right -> Coordinate`

Halo expansion on the positive (right) side of each direction.

#### `id_node -> int`

Linear index of this MPI process within the grid.

#### `num_node -> int`

Total number of MPI processes: `product(size_node)`.

#### `coor_node -> Coordinate`

This process's position in the MPI grid (e.g., `Coordinate([0, 0, 0, 0])`
for the first process).

#### `size_node -> Coordinate`

Dimensions of the MPI grid (e.g., `Coordinate([1, 1, 1, 2])` for 2 nodes
along the time direction).

#### `node_site -> Coordinate`

Local site dimensions per node. Same as `local_site`.

#### `is_only_local -> bool`

`True` if the geometry has no halo expansion (i.e., `expansion_left` and
`expansion_right` are both zero).

---

### Coordinate Mapping

#### `coordinate_g_from_l(xl: Coordinate) -> Coordinate`

Convert a local coordinate `xl` to a global coordinate `xg`:

```
xg = xl + coor_node * node_site
```

#### `coordinate_l_from_g(xg: Coordinate) -> Coordinate`

Convert a global coordinate `xg` to a local coordinate `xl`:

```
xl = xg - coor_node * node_site
```

#### `index_from_coordinate(xl: Coordinate) -> int`

Map a local coordinate to a linear site index. Returns `0 <= index <
local_volume()`.

#### `coordinate_from_index(index: int) -> Coordinate`

Map a linear site index back to a local coordinate. Inverse of
`index_from_coordinate`.

#### `index_from_g_coordinate(xg: Coordinate) -> int`

Map a global coordinate to a local linear site index. Combines
`coordinate_l_from_g` and `index_from_coordinate`.

#### `g_coordinate_from_index(index: int) -> Coordinate`

Map a local linear site index to a global coordinate. Combines
`coordinate_from_index` and `coordinate_g_from_l`.

---

### Site Queries

#### `is_local(xl: Coordinate) -> bool`

Return whether the local coordinate `xl` falls within the unexpanded local
volume (`0 <= xl[mu] < node_site[mu]` for each direction, or `size_node[mu]
== 1`). Also checks even/odd consistency if `eo != 0`.

#### `is_local_xg(xg: Coordinate) -> bool`

Return whether the global coordinate `xg` falls within this node's local
volume. Internally converts to a local coordinate and calls `is_local`.

#### `xg_arr() -> numpy.ndarray`

Return a NumPy array of shape `(local_volume, 4)` containing the global
coordinate of every local site. Requires `field_selection` module.

---

### Display and Pickle

#### `__str__` / `__repr__` / `show()`

Return a human-readable string. For a plain geometry:
```
Geometry([4, 4, 4, 8])
```
For a geometry with expansion or even/odd:
```
Geometry([4, 4, 4, 8], expan_left=[1, 1, 1, 1], expan_right=[1, 1, 1, 1], eo=0)
```

#### `__getstate__` / `__setstate__`

Pickle support. The state is node-dependent: it includes `size_node`,
`coor_node`, `node_site`, `expansion_left`, and `expansion_right`. Pickling
and unpickling on different MPI grids may produce different geometries.

---

## `geo_resize` Function

```python
geo_resize(geo, expansion_left=None, expansion_right=None) -> Geometry
```

Create a copy of `geo` with the specified halo expansion. This is the
primary way to add halo regions to a geometry for stencil operations
(e.g., computing gauge links that require neighboring sites).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `geo` | `Geometry` | — | Source geometry |
| `expansion_left` | `None`, `int`, `list`, or `Coordinate` | `None` (→ 0) | Halo expansion on the negative side |
| `expansion_right` | `None`, `int`, `list`, or `Coordinate` | `None` (→ 0) | Halo expansion on the positive side |

When an `int` is given, the same thickness is applied to all four directions.
When a `list` is given, it is converted to a `Coordinate`.

**Note:** When `size_node[mu] == 1` (i.e., the lattice is not partitioned in
direction `mu`), the expansion in that direction is forced to zero because no
halo exchange is needed. On a single-node layout (`size_node = [1,1,1,1]`),
all expansions are zero.

```python
geo_base = q.Geometry(q.Coordinate([4, 4, 4, 8]))
geo_expanded = q.geo_resize(geo_base, 1)                # expand both sides by 1
geo_asymmetric = q.geo_resize(geo_base, [1,1,1,0], [1,1,1,0])  # no time-direction expansion
```

---

## `geo_eo` Function

```python
geo_eo(geo, eo=0) -> Geometry
```

Create a copy of `geo` with the even/odd flag set to `eo`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `geo` | `Geometry` | — | Source geometry |
| `eo` | `int` | `0` | `0` = full lattice, `1` = odd sites, `2` = even sites |

```python
geo_base = q.Geometry(q.Coordinate([4, 4, 4, 8]))
geo_odd = q.geo_eo(geo_base, 1)
print(geo_odd.eo)              # 1
print(geo_odd.local_volume)    # 256 (half of 512)
```

---

## Examples

### Basic Geometry Construction

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

print(f"total_site: {geo.total_site}")         # [4, 4, 4, 8]
print(f"local_site: {geo.local_site}")         # depends on num_node
print(f"total_volume: {geo.total_volume}")     # 512
print(f"local_volume: {geo.local_volume}")     # 512 / num_node
print(f"size_node: {geo.size_node}")           # depends on MPI config
print(f"id_node: {geo.id_node}")               # 0 .. num_node-1

q.end_with_mpi()
```

### Coordinate Mapping

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Local coordinate to global
xl = q.Coordinate([1, 2, 3, 4])
xg = geo.coordinate_g_from_l(xl)
print(f"local {xl.to_list()} -> global {xg.to_list()}")

# Global coordinate to local
xl_back = geo.coordinate_l_from_g(xg)
assert xl == xl_back

# Coordinate <-> index
idx = geo.index_from_coordinate(xl)
xl_from_idx = geo.coordinate_from_index(idx)
assert xl == xl_from_idx

q.end_with_mpi()
```

### Geometry with Halo Expansion

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo_base = q.Geometry(total_site)
geo_exp = q.geo_resize(geo_base, 1)

print(f"local_volume: {geo_exp.local_volume}")                # 512
print(f"local_volume_expanded: {geo_exp.local_volume_expanded}")
print(f"expansion_left: {geo_exp.expansion_left.to_list()}")
print(f"expansion_right: {geo_exp.expansion_right.to_list()}")
print(f"is_only_local: {geo_exp.is_only_local}")

# On a single-node layout, expansion is forced to zero in all directions
# because size_node[mu] == 1 means no halo exchange is needed.
# Output:
#   expansion_left:  [0, 0, 0, 0]
#   expansion_right: [0, 0, 0, 0]
#   is_only_local:   True

# On a multi-node layout (e.g., size_node = [1, 1, 1, 2]):
#   expansion_left:  [1, 1, 1, 1]
#   expansion_right: [1, 1, 1, 1]
#   is_only_local:   False

q.end_with_mpi()
```

### Even/Odd Geometry

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

geo_odd = q.geo_eo(geo, 1)
geo_even = q.geo_eo(geo, 2)

print(f"full local_volume: {geo.local_volume}")       # 512
print(f"odd  local_volume: {geo_odd.local_volume}")   # 256
print(f"even local_volume: {geo_even.local_volume}")  # 256

q.end_with_mpi()
```

### Pickle Round-Trip

```python
import qlat as q
import pickle

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

data = pickle.dumps(geo)
geo2 = pickle.loads(data)
assert geo == geo2

# With expansion
geo_exp = q.geo_resize(geo, 1)
data2 = pickle.dumps(geo_exp)
geo_exp2 = pickle.loads(data2)
assert geo_exp == geo_exp2

q.end_with_mpi()
```
