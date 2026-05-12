# `qlat.mpi` — Low-Level MPI Utilities

Source: `qlat/qlat/mpi.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [MPI Communicator Lifecycle](#mpi-communicator-lifecycle)
   - [`begin`](#begin)
   - [`end`](#end)
3. [Node Layout Queries](#node-layout-queries)
   - [`get_size_node`](#get_size_node)
   - [`get_coor_node`](#get_coor_node)
4. [Global Sum Reductions](#global-sum-reductions)
   - [`glb_sum`](#glb_sum)
   - [Type-Specific Global Sums](#type-specific-global-sums)
5. [Broadcast Operations](#broadcast-operations)
   - [`bcast_long`, `bcast_double`, `bcast_complex`](#bcast_long-bcast_double-bcast_complex)
   - [`bcast_lat_data`](#bcast_lat_data)
6. [IO Shuffle Utilities](#io-shuffle-utilities)
   - [`get_id_node_list_for_shuffle`](#get_id_node_list_for_shuffle)
7. [Examples](#examples)

---

## Overview

`mpi` provides the low-level MPI communication primitives used throughout qlat:
begin/end of MPI communicators, node layout queries, global sum reductions, and
broadcasts. These are the building blocks that higher-level functions
(`begin_with_mpi`, `end_with_mpi`, `bcast_py`, etc.) rely on.

Most users should use the high-level wrappers `q.begin_with_mpi(size_node_list)`
and `q.end_with_mpi()` from `mpi_utils.py` rather than calling `begin`/`end`
directly. The `glb_sum` generic dispatcher is the recommended interface for
global sums.

---

## MPI Communicator Lifecycle

### `begin`

```python
begin(id_node, size_node, color=0)
```

Initialize an MPI communicator for this node. Low-level; prefer
`begin_with_mpi(size_node_list)` for user code.

| Parameter | Type | Description |
|---|---|---|
| `id_node` | `int` | Node index within the MPI grid |
| `size_node` | `Coordinate` | Dimensions of the MPI grid |
| `color` | `int` | MPI communicator color (default 0) |

### `end`

```python
end(is_preserving_cache=False)
```

Finalize the currently active MPI communicator. Low-level; prefer `end_with_mpi()`
for user code. Clears cached data unless `is_preserving_cache` is `True`.

---

## Node Layout Queries

### `get_size_node`

```python
get_size_node() -> Coordinate
```

Return the dimensions of the MPI grid as a `Coordinate`. For example, with 2
processes along the last dimension, `get_size_node()` returns
`Coordinate([1, 1, 1, 2])`.

### `get_coor_node`

```python
get_coor_node() -> Coordinate
```

Return this MPI process's position within the grid as a `Coordinate`. For
example, process 0 returns `Coordinate([0, 0, 0, 0])` and process 1 returns
`Coordinate([0, 0, 0, 1])` in a `1×1×1×2` grid.

---

## Global Sum Reductions

### `glb_sum`

```python
glb_sum(x) -> same_type_as_x
```

Compute the global MPI sum of `x` across all processes. The input is **not**
modified; the result is returned as a new value.

Dispatches by type:

| `type(x)` | Behavior |
|---|---|
| `float` | Calls `glb_sum_double` |
| `complex` | Calls `glb_sum_complex` |
| `int` / `np.int64` | Calls `glb_sum_long` |
| `np.ndarray` | Calls `glb_sum_np` |
| `list` | Recursively sums each element |
| `tuple` | Recursively sums each element |
| `LatData` | Calls `ld.glb_sum()` |

### Type-Specific Global Sums

#### `glb_sum_long(x) -> int`

Global sum of a long integer. Modifies `x` in place (for C++ interop).
Recommend using `glb_sum(x)` instead for a non-mutating interface.

#### `glb_sum_double(x) -> float`

Global sum of a double. Same caveat about in-place modification.

#### `glb_sum_complex(x) -> complex`

Global sum of a complex number. Same caveat about in-place modification.

#### `glb_sum_np(x) -> np.ndarray`

Global sum of a NumPy array of any shape and dtype (`float64`, `int64`, or
`complex128`). Returns a new array; `x` is **not** modified.

#### `glb_sum_lat_data(ld) -> LatData`

Global sum of `LatData`. Returns a copy with the sum applied.

---

## Broadcast Operations

### `bcast_long`, `bcast_double`, `bcast_complex`

```python
bcast_long(x, root=0) -> int
bcast_double(x, root=0) -> float
bcast_complex(x, root=0) -> complex
```

Broadcast a scalar value from the `root` process to all processes. The input
`x` is modified in place on non-root processes and returned.

### `bcast_lat_data`

```python
bcast_lat_data(ld, root=0) -> LatData
```

Broadcast a `LatData` object from the `root` process. Returns a copy of `ld`
with the broadcast applied (the original `ld` is not modified).

---

## IO Shuffle Utilities

### `get_id_node_list_for_shuffle`

```python
get_id_node_list_for_shuffle() -> list[int]
```

Return the list of node indices used for distributed IO shuffling. The returned
list maps `id_node_in_shuffle → id_node`, i.e.,
`list[id_node_in_shuffle] = id_node`.

---

## Examples

### Global Sum of Scalars and Arrays

```python
import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
]

q.begin_with_mpi(size_node_list)

id_node = q.get_id_node()

# Scalar global sums
a = 12 + id_node
total = q.glb_sum(a)         # e.g., 25 with 2 nodes

b = 12.4 + id_node
total_b = q.glb_sum(b)       # e.g., 25.4 with 2 nodes

# NumPy array global sum
arr = np.arange(3.0) + 1.0 + id_node
total_arr = q.glb_sum(arr)   # element-wise sum across nodes

# Nested list global sum
lst = [1.0 + id_node, 2.0 + 1.0j * id_node]
total_lst = q.glb_sum(lst)   # recursive element-wise sum

q.end_with_mpi()
```

### Broadcast from Root

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
]

q.begin_with_mpi(size_node_list)

id_node = q.get_id_node()

# Root prepares the value; all nodes receive it
if id_node == 0:
    val = 42
else:
    val = 0
val = q.bcast_long(val)
# val == 42 on all nodes

q.end_with_mpi()
```

### Node Layout Inspection

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
]

q.begin_with_mpi(size_node_list)

print(q.get_id_node())       # 0 or 1
print(q.get_num_node())      # 2
print(q.get_coor_node())     # Coordinate([0, 0, 0, 0]) or Coordinate([0, 0, 0, 1])
print(q.get_size_node())     # Coordinate([1, 1, 1, 2])

q.end_with_mpi()
```
