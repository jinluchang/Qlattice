# `qlat.psel_split` — Point-Selection Splitting with Maximum Separation

Source: `qlat/qlat/psel_split.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Key Concepts](#key-concepts)
3. [`PointsDistanceTree` Class](#pointsdistancetree-class)
4. [`PointsDistanceSet` Class](#pointsdistanceset-class)
5. [Module-Level Functions](#module-level-functions)
   - [Closest-Point Search](#closest-point-search)
   - [Splitting Functions](#splitting-functions)
6. [Examples](#examples)

---

## Overview

`psel_split` partitions a `PointsSelection` (a set of lattice points) into
two or more subsets while **maximising the minimum separation** between
points that share the same subset. This is useful when independent
measurements must be taken at well-separated lattice sites to reduce
autocorrelation.

The module builds a spatial tree (`PointsDistanceTree`) that stores points
grouped by distance ranges, enabling efficient nearest-neighbour queries.
Two splitting strategies are provided:

- **closest** — greedily assigns each point to the subset whose nearest
  already-assigned neighbour is farthest away.
- **ranking** — assigns based on a configurable ranking function over the
  *n* closest neighbours (default, more robust).

---

## Key Concepts

- **Periodic distance** — all distances are computed with periodic boundary
  conditions via `smod_coordinate`, so the lattice wraps around.
- **Spatial tree** — points are stored in a hierarchy of distance ranges.
  Each node holds a representative point and a `range_sqr` value; children
  lie within that squared-distance radius. This allows pruning entire
  subtrees during nearest-neighbour searches.
- **MPI parallelism** — the heavy `find_all_closest_point_list` and
  `find_all_closest_n_point_list` routines distribute work across MPI ranks
  via `get_mpi_chunk` and `parallel_map`.

---

## `PointsDistanceTree` Class

A recursive spatial tree that stores lattice points for efficient
distance queries.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `point` | `Coordinate` or `None` | Representative point of this node |
| `range_sqr` | `int` or `None` | Squared-distance radius (leaf if `None`) |
| `tree_list` | `list` or `None` | Child subtrees |

### Class Methods

#### `mk_with_point(xg, range_sqr=None) -> PointsDistanceTree`

Create a leaf node for a single point (when `range_sqr` is `None` or `0`),
or a tree node enclosing `xg` with at least the given `range_sqr`.

#### `mk_with_tree(t0) -> PointsDistanceTree`

Create a new tree node one level above `t0`, doubling the range.

### Instance Methods

#### `add(xg, total_site) -> PointsDistanceTree`

Return a (possibly new) tree with point `xg` inserted. May modify the
current tree in place.

#### `try_add(xg, total_site) -> bool`

Attempt to add `xg` in place. Returns `True` on success.

#### `check(total_site)`

Assert internal consistency of the tree structure.

#### `count() -> int`

Return the number of points stored in the tree.

#### `list() -> list`

Return a nested list representation of the tree.

#### `find_closest_point_list(xg, total_site) -> (int, list)`

Find all points in the tree at the minimum squared distance from `xg`.
Returns `(mini_dis_sqr, point_list)`. The point `xg` itself is excluded.

#### `find_closest_n_point_list(xg, n, total_site) -> list`

Find the `n` closest points to `xg`. Returns a list of `(dis_sqr, xg)`
tuples sorted by distance. The point `xg` itself is excluded.

---

## `PointsDistanceSet` Class

A convenience wrapper that builds a `PointsDistanceTree` from a
`PointsSelection`.

### Constructor

```python
PointsDistanceSet(psel: PointsSelection = None)
```

### Methods

#### `set_psel(psel)`

Build the internal tree from a `PointsSelection`.

#### `add(xg)`

Insert an additional point.

#### `check()`

Assert internal consistency.

#### `count() -> int`

Return the number of stored points.

#### `list() -> list`

Return a nested list representation.

#### `find_closest_point_list(xg) -> (int, list)`

Find all points at the minimum squared distance from `xg` (excludes `xg`).

#### `find_closest_n_point_list(xg, n) -> list`

Find the `n` closest points to `xg`.

---

## Module-Level Functions

### Closest-Point Search

#### `find_all_closest_point_list(psel, rs=None, is_parallel=True) -> list`

For every point in `psel`, find its closest neighbour(s). Returns a list of
`(mini_dis_sqr, (xg, [xg1, ...]))` tuples, sorted by distance and randomly
permuted to break ties.

#### `find_all_closest_n_point_list(psel, n, ranking_func=None, rs=None) -> list`

For every point in `psel`, find the `n` closest neighbours and rank them.
Returns a list of `(ranking, (xg, [(dis_sqr, xg1), ...]))` tuples, sorted
by ranking.

#### `find_all_closest_n_point_list_ranking_func_default(dis_sqr_list) -> float`

Default ranking function: computes a combined score from the list of squared
distances using a power-law weighting.

#### `find_closest_dis_sqr_for_psel_list(psel_list, is_parallel=True) -> list`

For each `PointsSelection` in `psel_list`, find the minimum closest-point
squared distance. Useful for verifying separation after splitting.

### Splitting Functions

#### `psel_split_that_increase_separation(psel, mode=None, rs=None) -> (PointsSelection, PointsSelection)`

Split `psel` into two subsets. `mode` is `"ranking"` (default) or `"closest"`.

#### `psel_split_that_increase_separation_closest(psel, rs=None) -> (PointsSelection, PointsSelection)`

Split using the closest-neighbour greedy strategy.

#### `psel_split_that_increase_separation_ranking(psel, n, ranking_func=None, rs=None) -> (PointsSelection, PointsSelection)`

Split using the ranking-based strategy over the `n` closest neighbours.

#### `psel_split_n_that_increase_separation(psel, num_piece, rs=None) -> list`

Recursively split `psel` into `num_piece` subsets. Returns a list of
`PointsSelection` objects of length `num_piece`.

---

## Examples

### Splitting a Point Selection in Two

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([8, 8, 8, 16])
xg_list = [
    q.Coordinate([0, 0, 0, 0]),
    q.Coordinate([1, 0, 0, 0]),
    q.Coordinate([4, 4, 4, 8]),
    q.Coordinate([4, 4, 4, 9]),
]
psel = q.PointsSelection(total_site, xg_list)

from qlat.psel_split import psel_split_that_increase_separation
psel1, psel2 = psel_split_that_increase_separation(psel)
print(f"psel1 has {len(psel1)} points, psel2 has {len(psel2)} points")

q.end_with_mpi()
```

### Splitting into N Pieces

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([8, 8, 8, 16])
xg_list = [q.Coordinate([i, 0, 0, 0]) for i in range(8)]
psel = q.PointsSelection(total_site, xg_list)

from qlat.psel_split import psel_split_n_that_increase_separation
psel_list = psel_split_n_that_increase_separation(psel, num_piece=4)
for i, p in enumerate(psel_list):
    print(f"subset {i}: {len(p)} points")

q.end_with_mpi()
```

### Querying Closest Points

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([8, 8, 8, 16])
xg_list = [
    q.Coordinate([0, 0, 0, 0]),
    q.Coordinate([2, 0, 0, 0]),
    q.Coordinate([0, 3, 0, 0]),
]
psel = q.PointsSelection(total_site, xg_list)

from qlat.psel_split import PointsDistanceSet
pds = PointsDistanceSet(psel)
query = q.Coordinate([0, 0, 0, 0])
dis_sqr, closest = pds.find_closest_point_list(query)
print(f"closest distance squared: {dis_sqr}")

q.end_with_mpi()
```
