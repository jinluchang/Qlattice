# `qlat.mpi_utils` — MPI Initialization and Data Distribution

Source: `qlat/qlat/mpi_utils.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [MPI Lifecycle Functions](#mpi-lifecycle-functions)
3. [Communication and Diagnostics](#communication-and-diagnostics)
4. [Data Distribution](#data-distribution)
5. [Default Node Grid](#default-node-grid)
6. [Examples](#examples)

---

## Overview

`mpi_utils` handles MPI setup and teardown for qlat. When qlat runs in
parallel it must know the 4-D Cartesian grid of MPI ranks (`size_node`);
this module maps the world communicator to such a grid automatically.

The key entry points are:

- `begin_with_mpi` — initialize MPI and the qlat C++ runtime.
- `end_with_mpi` — finalize both.
- `get_mpi_chunk` — split a list across MPI ranks for distributed work.

The module also stores a module-level `comm` object (an `mpi4py`
communicator) that other qlat modules use for collective operations.

---

## MPI Lifecycle Functions

### `begin_with_mpi(size_node_list=None)`

Initialize MPI and the qlat C++ runtime.

1. Calls `MPI.COMM_WORLD` to obtain the world communicator.
2. Looks up the number of MPI ranks (`comm.size`) and tries to match it
   against `size_node_list` (a list of `[n0, n1, n2, n3]` or `Coordinate`
   objects). The user-supplied list is prepended to
   `default_size_node_list`.
3. Calls `q.begin(id_node, size_node)` to initialize the C++ layer.

Raises `Exception` if no matching grid is found for the current `num_node`.

```python
begin_with_mpi(size_node_list=None)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `size_node_list` | `list` or `None` | `None` | Candidate 4-D node grids to try before the built-in defaults |

### `end_with_mpi(is_preserving_cache=False)`

Finalize the qlat C++ runtime (`q.end`) and call `MPI.Finalize()`.

```python
end_with_mpi(is_preserving_cache=False)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `is_preserving_cache` | `bool` | `False` | If `True`, preserve the qlat cache across calls |

---

## Communication and Diagnostics

### `get_comm() -> mpi4py.MPI.Intracomm`

Return the module-level MPI communicator set by `begin_with_mpi`.

### `set_comm(x)`

Set the module-level MPI communicator. Rarely called directly.

### `show_machine()`

Print a diagnostic line showing the current rank ID, total ranks, and
Cartesian coordinates of this node. Decorated with `@timer_verbose`.

---

## Data Distribution

### `get_mpi_chunk(total_list, *, rng_state=None) -> list`

Distribute `total_list` across MPI ranks. Each rank receives one chunk of
the list. The assignment is deterministic when `rng_state` is `None`; pass
an `RngState` for reproducible shuffling before chunking.

```python
get_mpi_chunk(total_list, *, rng_state=None)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `total_list` | `list` | — | The full list to distribute |
| `rng_state` | `RngState` or `None` | `None` | Optional RNG for shuffling (must be identical on all ranks) |

Returns the chunk assigned to the current rank (may be `[]` if there are
more ranks than items).

Internally uses `qlat_utils.get_chunk_list` with `chunk_number=get_num_node()`.

---

## Default Node Grid

`default_size_node_list` contains 27 pre-defined 4-D grids for common
node counts from 1 to 512. The list includes:

| Nodes | Grid (n0 x n1 x n2 x n3) |
|-------|---------------------------|
| 1     | 1 x 1 x 1 x 1            |
| 2     | 1 x 1 x 1 x 2            |
| 4     | 1 x 1 x 2 x 2            |
| 8     | 1 x 2 x 2 x 2            |
| 16    | 2 x 2 x 2 x 2            |
| 32    | 2 x 2 x 2 x 4            |
| 64    | 2 x 2 x 4 x 4            |
| 128   | 2 x 4 x 4 x 4            |
| 256   | 4 x 4 x 4 x 4            |
| 512   | 4 x 4 x 4 x 8            |
| ...   | (see source for full list) |

Grids with a factor of 3 in the last dimension (3, 6, 12) are also
included for anisotropic lattice geometries.

---

## Examples

### Basic MPI Setup and Teardown

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

q.mpi_utils.show_machine()

q.end_with_mpi()
```

### Distributing Work Across Ranks

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

full_list = list(range(100))
my_chunk = q.mpi_utils.get_mpi_chunk(full_list)
print(f"Rank {q.get_id_node()} processing {len(my_chunk)} items")

q.end_with_mpi()
```

### Custom Node Grid

```python
import qlat as q

# Use a 2x2x2x4 grid for 32 ranks
q.begin_with_mpi([[2, 2, 2, 4]])

print(f"size_node = {q.get_size_node()}")

q.end_with_mpi()
```
