# `qlat.field_selection` — Lattice Site Selection and Shuffle Plans

Source: `qlat/qlat/field_selection.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`PointsSelection` Class](#pointsselection-class)
   - [Constructors](#constructors)
   - [Properties](#properties)
   - [Indexing and Iteration](#indexing-and-iteration)
   - [Set Operations](#set-operations)
   - [Serialization and I/O](#serialization-and-io)
   - [MPI Broadcast](#mpi-broadcast)
   - [Hashing and Pickle](#hashing-and-pickle)
3. [`FieldSelection` Class](#fieldselection-class)
   - [Constructors](#fieldselection-constructors)
   - [Properties](#fieldselection-properties)
   - [Selection Modification](#selection-modification)
   - [Set Operations](#fieldselection-set-operations)
   - [Indexing and Iteration](#fieldselection-indexing-and-iteration)
   - [I/O and Pickle](#fieldselection-io-and-pickle)
4. [`SelectedShufflePlan` Class](#selectedshuffleplan-class)
   - [Constructors](#selectedshuffleplan-constructors)
   - [Properties](#selectedshuffleplan-properties)
   - [Shuffle Methods](#shuffle-methods)
5. [Module-Level Functions](#module-level-functions)
6. [Examples](#examples)

---

## Overview

`field_selection` provides the mechanism for working with subsets of lattice
sites in a distributed (MPI) environment. Two complementary selection classes
are defined:

- **`PointsSelection`** — a list of global coordinates (`xg`) representing
  selected sites. Supports multiple distribution types: Global (`"g"`),
  Full (`"f"`), Local (`"l"`), Random (`"r"`), and Other (`"o"`).
- **`FieldSelection`** — a per-site rank array stored as a `Field`. A site is
  selected if its rank is `>= 0`. Provides efficient lookup from global
  coordinate to local index.

**`SelectedShufflePlan`** builds communication plans for redistributing
selected points between different distribution types (e.g., Local to Global,
Local to Random, or time-slice shuffles). It is used internally by
`SelectedField` and `SelectedPoints` operations.

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
psel = q.PointsSelection(total_site)
psel.set_rand(total_site, 10, q.RngState("seed"))
print(psel.n_points)         # 10
print(psel.points_dist_type) # "g"

fsel = q.FieldSelection(psel, q.Geometry(total_site))
print(fsel.n_elems)          # 10

q.end_with_mpi()
```

---

## `PointsSelection` Class

A selection of lattice sites represented as a list of 4-D global coordinates.

### Constructors

### `PointsSelection()`

Create an empty, uninitialized selection.

### `PointsSelection(total_site: Coordinate)`

Create a zero-length selection with the given total lattice dimensions.

| Parameter | Type | Description |
|---|---|---|
| `total_site` | `Coordinate` | Total lattice dimensions |

### `PointsSelection(total_site, xg_arr)`

Create from an explicit list of coordinates.

| Parameter | Type | Description |
|---|---|---|
| `total_site` | `Coordinate` | Total lattice dimensions |
| `xg_arr` | `int`, `Coordinate`, `numpy.ndarray`, or `list` | If `int`, treated as `n_points`; if `Coordinate`, a single point; if `ndarray`, shape `(n_points, 4)`; if `list` of `Coordinate` |

### `PointsSelection(total_site, xg_arr, points_dist_type)`

Same as above but also set the distribution type.

| Parameter | Type | Description |
|---|---|---|
| `points_dist_type` | `str` | One of `"g"`, `"f"`, `"l"`, `"r"`, `"o"` |

### `PointsSelection(geo: Geometry)`

Create a full selection (all sites) with `points_dist_type = "f"`.

### `PointsSelection(fsel: FieldSelection)`

Create a local selection from a `FieldSelection` (`points_dist_type = "l"`).

### `PointsSelection(fsel, ssp)`

Create by shuffling a `FieldSelection` according to a `SelectedShufflePlan`.

### `PointsSelection(psel: PointsSelection)`

Copy constructor.

### `PointsSelection(psel, ssp, is_reverse=False)`

Create by shuffling an existing `PointsSelection` according to a
`SelectedShufflePlan`.

---

### Properties

### `points_dist_type -> str`

Get or set the distribution type. Values:
- `"g"` — Global: all points known on every node.
- `"f"` — Full: every lattice site selected.
- `"l"` — Local: points distributed across nodes (each node holds its local subset).
- `"r"` — Random: points randomly distributed across nodes.
- `"o"` — Other.

### `total_site -> Coordinate`

Get or set the total lattice dimensions.

### `geo -> Geometry`

A `Geometry` constructed from `total_site` with no halo expansion.

### `n_points -> int`

Number of selected points on this node.

### `xg_arr -> numpy.ndarray`

The global coordinates of all selected points as a `(n_points, 4)` array of
`int32`. Setting this property accepts an `int` (number of points), a single
`Coordinate`, an `ndarray`, or a `list` of coordinates.

---

### Indexing and Iteration

### `__getitem__(idx)` / `__setitem__(idx, val)`

Index into the underlying `xg_arr` via `numpy.asarray`.

### `__iter__()`

Iterate over selected points, yielding `Coordinate` objects.

### `__len__()`

Return `n_points`.

### `coordinate_from_idx(idx: int) -> Coordinate`

Return the global coordinate at linear index `idx`.

---

### Set Operations

### `intersect(fsel: FieldSelection) -> PointsSelection`

Return a new `PointsSelection` containing only points that are also selected
by `fsel`.

### `is_containing(sel_small) -> bool`

Return `True` if every point in `sel_small` (a `PointsSelection` or
`FieldSelection`) is contained in `self`.

### `is_containing_psel(psel_small: PointsSelection) -> bool`

Check containment for a `PointsSelection`.

### `is_containing_fsel(fsel_small: FieldSelection) -> bool`

Check containment for a `FieldSelection`.

---

### Serialization and I/O

### `save(path: str, *, is_sync_node=True)`

Save to file. When `is_sync_node=True`, writes via synchronized I/O
(`save_points_selection_info`); otherwise writes per-node
(`save_points_selection`).

### `load(path: str, geo=None, *, is_sync_node=True)`

Load from file. For `.lati` files, `geo` is optional (total_site is stored
in the file). For other formats, `geo` is required.

### `save_str() -> bytes`

Serialize to a byte string (only node 0 returns data).

### `load_str(content: bytes)`

Deserialize from a byte string (only node 0 needs the data; result is
broadcast).

### `to_lat_data() -> LatDataInt`

Convert to `LatDataInt` format.

### `from_lat_data(ld: LatDataInt)`

Load from `LatDataInt`.

---

### MPI Broadcast

### `bcast(root=0) -> PointsSelection`

Broadcast the selection from `root` to all nodes using native MPI broadcast.

### `bcast_via_ld(root=0) -> PointsSelection`

Broadcast via `LatDataInt` serialization (alternative implementation).

---

### Hashing and Pickle

### `hash_sha256() -> str`

Return a SHA-256 hash that is consistent across all nodes. For `"g"`
distribution, returns the local hash if all nodes agree; otherwise hashes
the gathered data from all nodes.

### `__getstate__` / `__setstate__`

Pickle support. Only works correctly on a single node or when all nodes hold
the same data.

---

## `FieldSelection` Class

A site-level selection stored as a rank field. Each local site has a rank
value: `rank >= 0` means selected, `rank == -1` means not selected.

### <a id="fieldselection-constructors"></a> Constructors

### `FieldSelection()`

Create an empty, uninitialized selection.

### `FieldSelection(geo: Geometry, rank=-1)`

Create a uniform selection. `rank=0` selects all sites; `rank=-1` (default)
selects no sites.

### `FieldSelection(psel: PointsSelection, geo=None)`

Create from a `PointsSelection`. Requires `psel.points_dist_type in ["l", "f", "g"]`.

---

### <a id="fieldselection-properties"></a> Properties

### `geo -> Geometry`

The geometry of the underlying rank field.

### `total_site -> Coordinate`

The total lattice dimensions.

### `n_elems -> int`

Number of selected elements on this node.

---

### Selection Modification

### `update()`

Rebuild internal indices from the `f_rank` field. Must be called after
modifying `f_rank` directly (e.g., via buffer view).

### `set_empty(geo: Geometry)`

Set an empty selection (all ranks = -1) with the given geometry.

### `set_uniform(geo: Geometry, val=0)`

Set a uniform selection. `val=0` selects all sites; `val=-1` deselects all.

### `set_rand(total_site: Coordinate, n_per_tslice: int, rs: RngState)`

Randomly select `n_per_tslice` sites per time slice.

### `set_rand_psel(total_site, n_per_tslice, rs, psel=None)`

Same as `set_rand`, then additionally include all points from `psel`.

### `add_psel(psel: PointsSelection, rank_psel=...)`

Add points from a `PointsSelection` to the selection. Points already selected
with a lower rank keep their existing rank.

### `add_fsel(fsel: FieldSelection)`

Add points from another `FieldSelection`. Points already selected with a
lower rank keep their existing rank.

---

### <a id="fieldselection-set-operations"></a> Set Operations

### `intersect_with(fsel: FieldSelection)`

Modify `self` to contain only sites also selected by `fsel`. More efficient
if `self` is smaller.

### `intersect(fsel: FieldSelection) -> FieldSelection`

Return a new `FieldSelection` containing only sites selected by both.
Does not modify `self`.

### `is_containing(sel_small) -> bool`

Return `True` if every point in `sel_small` is contained in `self`.

### `is_containing_psel(psel_small) -> bool` / `is_containing_fsel(fsel_small) -> bool`

Specialized containment checks.

### `to_psel() -> PointsSelection`

Convert to a Global-distribution `PointsSelection`.

### `to_psel_local() -> PointsSelection`

Convert to a Local-distribution `PointsSelection`.

---

### <a id="fieldselection-indexing-and-iteration"></a> Indexing and Iteration

### `__getitem__(idx)` / `__setitem__(idx, val)`

Index into the rank array via `numpy.asarray`. Modifying values requires
calling `update()` afterwards.

### `__iter__()`

Iterate over selected sites, yielding `Coordinate` objects.

### `__len__()`

Return `n_elems`.

### `idx_from_coordinate(xg: Coordinate) -> int`

Look up the local index for a global coordinate. Returns the index in the
selection (not the field index).

### `coordinate_from_idx(idx: int) -> Coordinate`

Return the global coordinate for the selected element at index `idx`.

---

### <a id="fieldselection-io-and-pickle"></a> I/O and Pickle

### `save(path: str) -> int`

Write the selection to a file. Returns total bytes written.

### `load(path: str) -> int`

Read the selection from a file. Returns total bytes read.

### `__getstate__` / `__setstate__`

Pickle support. Only works on a single node.

---

## `SelectedShufflePlan` Class

Builds an MPI communication plan for redistributing `SelectedPointsChar`
data between different point distribution types. The plan is reusable and
can shuffle forward or in reverse.

### <a id="selectedshuffleplan-constructors"></a> Constructors

### `SelectedShufflePlan()`

Create an empty (identity) plan.

### `SelectedShufflePlan("l_from_g", psel_src, root)`

Shuffle from Global to Local distribution. `root` is the node that holds
the full global point list.

### `SelectedShufflePlan("l_from_g", psel_src_list, root_list)`

Batch version: multiple point selections with corresponding root nodes.

### `SelectedShufflePlan("g_from_l", psel_src, root, geo)`

Shuffle from Local to Global distribution.

### `SelectedShufflePlan("g_from_l", psel_src_list, root_list, geo_src_list)`

Batch version.

### `SelectedShufflePlan("r_from_l", psel_src, geo, rs)`

Shuffle from Local to Random distribution.

### `SelectedShufflePlan("r_from_l", psel_src_list, geo_src_list, rs)`

Batch version.

### `SelectedShufflePlan("dist_r_from_l", psel_src, geo, rs, id_node_list)`

Shuffle from Local to Random, restricted to nodes in `id_node_list`.

### `SelectedShufflePlan("t_slice_from_l", psel_src_list, geo_src_list)`

Shuffle from Local to Local by time slice (transpose data across time
slices).

### `SelectedShufflePlan("dist_t_slice_from_l", psel_src, geo, num_field)`

Distributed time-slice shuffle. Matches `"t_slice_from_l"` for use cases
like spatial smearing where propagators and gauge fields need compatible
shuffles.

---

### <a id="selectedshuffleplan-properties"></a> Properties

### `points_dist_type_send -> str`

Distribution type of the source data.

### `points_dist_type_recv -> str`

Distribution type of the destination data.

### `num_selected_points_send -> int`

Number of points on the send side.

### `num_selected_points_recv -> int`

Number of points on the receive side.

### `psel_send_list -> list` / `psel_recv_list -> list`

Lists of `PointsSelection` objects for send and receive sides.

### `geo_send_list -> list` / `geo_recv_list -> list`

Lists of `Geometry` objects for send and receive sides (computed lazily).

### `fsel_send_list -> list` / `fsel_recv_list -> list`

Lists of `FieldSelection` objects derived from the point selections
(computed lazily; only populated for Local/Full distribution types).

---

### Shuffle Methods

### `shuffle(sp_src: SelectedPointsChar, *, is_reverse=False) -> SelectedPointsChar`

Shuffle a single `SelectedPointsChar` using this plan. If `is_reverse=True`,
shuffles in the reverse direction.

### `shuffle_list(sp_src_list, *, is_reverse=False) -> list`

Shuffle a list of `SelectedPointsChar` objects.

### `shuffle_sp(cls, src, *, is_reverse=False)`

Shuffle a single field object (`SelectedPointsBase`, `SelectedFieldBase`, or
`FieldBase`). Converts to `SelectedPointsChar` internally, shuffles, and
converts back to type `cls`.

### `shuffle_sp_list(cls, src_list, *, is_reverse=False)`

Shuffle a list of field objects.

---

## Module-Level Functions

### `mk_xg_field(geo: Geometry) -> FieldInt`

Create a `FieldInt` containing the global coordinate index for each local
site.

### `get_psel_single(total_site, xg=None) -> PointsSelection`

Return a cached single-point selection. If `xg` is `None`, uses
`[-1, -1, -1, -1]`.

### `get_psel_tslice(total_site, *, t_dir=3) -> PointsSelection`

Return a cached time-slice selection. For `t_dir=3`, selects all spatial
sites at each time value. For `t_dir=2`, selects along the z-direction.

### `is_matching_fsel(fsel1, fsel2) -> bool`

Check if two `FieldSelection` objects have matching selection patterns.

---

## Examples

### Creating a PointsSelection

```python
import qlat as q
import numpy as np

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])

# From a numpy array of coordinates
xg_arr = np.array([[0, 0, 0, t] for t in range(8)], dtype=np.int32)
psel = q.PointsSelection(total_site, xg_arr)
print(f"points_dist_type: {psel.points_dist_type}")  # "g"
print(f"n_points: {psel.n_points}")                  # 8

# From a single coordinate
psel_single = q.PointsSelection(total_site, q.Coordinate([1, 2, 3, 4]))
print(f"n_points: {psel_single.n_points}")            # 1

q.end_with_mpi()
```

### Creating a FieldSelection

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Select all sites
fsel_all = q.FieldSelection(geo, 0)
print(f"n_elems: {fsel_all.n_elems}")   # 512

# Select no sites
fsel_none = q.FieldSelection(geo)
print(f"n_elems: {fsel_none.n_elems}")  # 0

# Random selection
rs = q.RngState("test")
fsel_rand = q.FieldSelection()
fsel_rand.set_rand(total_site, 4, rs)
print(f"n_elems: {fsel_rand.n_elems}")  # 4 * 8 = 32

q.end_with_mpi()
```

### Using PointsSelection with a Geometry

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Full selection (every site)
psel_full = q.PointsSelection(geo)
print(f"points_dist_type: {psel_full.points_dist_type}")  # "f"
print(f"n_points: {psel_full.n_points}")                  # 512

q.end_with_mpi()
```

### Time-Slice Selection

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])

# Get a selection of all spatial sites per time slice
psel_tslice = q.get_psel_tslice(total_site)
print(f"n_points: {psel_tslice.n_points}")  # 8 (one entry per time slice)

q.end_with_mpi()
```

### Intersecting Selections

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Create a random field selection
rs = q.RngState("seed")
fsel = q.FieldSelection()
fsel.set_rand(total_site, 2, rs)

# Create a global point selection
psel = q.PointsSelection(total_site)
psel.set_rand(total_site, 20, q.RngState("seed2"))

# Intersect: keep only points in psel that are also in fsel
psel_intersect = psel.intersect(fsel)
print(f"psel n_points: {psel.n_points}")
print(f"fsel n_elems: {fsel.n_elems}")
print(f"psel_intersect n_points: {psel_intersect.n_points}")

q.end_with_mpi()
```

### Building a Shuffle Plan

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Create a local point selection from a random field selection
rs = q.RngState("test")
fsel = q.FieldSelection()
fsel.set_rand(total_site, 2, rs)
psel_local = q.PointsSelection(fsel)

# Build a shuffle plan: local -> global
ssp = q.SelectedShufflePlan("g_from_l", psel_local, 0, geo)
print(f"send dist type: {ssp.points_dist_type_send}")
print(f"recv dist type: {ssp.points_dist_type_recv}")

q.end_with_mpi()
```

### Saving and Loading a PointsSelection

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
rs = q.RngState("seed")
psel = q.PointsSelection(total_site)
psel.set_rand(total_site, 4, rs)

# Save to file
psel.save("/tmp/test_psel.lati")

# Load from file
psel_loaded = q.PointsSelection()
psel_loaded.load("/tmp/test_psel.lati", q.Geometry(total_site))
assert psel == psel_loaded

q.end_with_mpi()
```
