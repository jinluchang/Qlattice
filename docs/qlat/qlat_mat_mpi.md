# `qlat.mat_mpi` — Distributed NumPy Arrays over MPI

Source: `qlat/qlat/mat_mpi.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Module-Level Configuration](#module-level-configuration)
   - [`set_mpi_comm`](#set_mpi_comm)
   - [`get_mpi_comm`](#get_mpi_comm)
   - [`bcast_py`](#bcast_py)
3. [The `DistArray` Class](#the-distarray-class)
   - [Constructor](#constructor)
   - [Arithmetic Operators](#arithmetic-operators)
   - [`sum`](#distarraysum)
   - [`transpose` / `transpose2d`](#distarraytranspose--transpose2d)
   - [`conj`](#distarrayconj)
4. [Distributed Linear Algebra](#distributed-linear-algebra)
   - [`d_matmul`](#d_matmul)
   - [`d_trace`](#d_trace)
5. [Scatter and Gather](#scatter-and-gather)
   - [`scatter_arr`](#scatter_arr)
   - [`gather_arr`](#gather_arr)
   - [`all_gather_arr`](#all_gather_arr)
6. [Examples](#examples)

---

## Overview

The `qlat.mat_mpi` module provides a lightweight distributed-array layer on
top of NumPy and mpi4py. The central abstraction is `DistArray`, which
partitions the **first dimension** of a NumPy array evenly across MPI ranks
(with zero-padding when the size does not divide evenly).

Features:

- Element-wise arithmetic (`+`, `-`, `*`, `/`) between `DistArray` objects
  and between `DistArray` and `np.ndarray`.
- Distributed `sum`, `transpose`, `conj`, `matmul`, and `trace` operations.
- `scatter_arr` / `gather_arr` / `all_gather_arr` for moving data between
  a single `np.ndarray` and a `DistArray`.
- A `use_reference_implementation` flag for switching between optimized
  and reference (gather-to-root) implementations for debugging.

---

## Module-Level Configuration

### `set_mpi_comm`

```python
set_mpi_comm(comm) -> None
```

Set the default MPI communicator used when no `comm` argument is passed.
Must be called before any `DistArray` creation if a communicator other than
`MPI.COMM_WORLD` is desired.

---

### `get_mpi_comm`

```python
get_mpi_comm() -> mpi4py.MPI.Intracomm
```

Return the default MPI communicator. Falls back to `MPI.COMM_WORLD` if
`set_mpi_comm` has not been called.

---

### `bcast_py`

```python
bcast_py(x, root=0, comm=None) -> Any
```

Broadcast a Python object `x` from `root` to all ranks. A thin wrapper
around `comm.bcast`.

---

## The `DistArray` Class

### Constructor

```python
DistArray(*, comm=None)
```

Create an empty distributed array. The array is initialized to a single
`float64` zero. After construction, set `self.x` (local NumPy array) and
`self.n` (total first-dimension size) before use.

| Attribute | Type | Description |
|---|---|---|
| `n` | `int` | Total size of the distributed (first) dimension |
| `x` | `np.ndarray` | Local portion of the array (padded with zeros if needed) |
| `comm` | `MPI.Intracomm` | MPI communicator |

---

### Arithmetic Operators

All standard arithmetic operators are supported between two `DistArray`
objects (requiring the same communicator, total size, and number of
dimensions) and between `DistArray` and scalars / `np.ndarray`.

| Operator | Method | Description |
|---|---|---|
| `+` | `__add__`, `__radd__` | Element-wise addition |
| `-` | `__sub__`, `__rsub__` | Element-wise subtraction |
| `*` | `__mul__`, `__rmul__` | Element-wise multiplication |
| `/` | `__truediv__`, `__rtruediv__` | Element-wise division |
| `@` | `__matmul__` | Distributed matrix-vector product (delegates to `d_matmul`) |

All arithmetic operations return a new `DistArray`.

---

### `DistArray.sum`

```python
sum(axis=None, *, keepdims=False) -> np.ndarray
```

Compute the sum across all ranks. Collective operation. The result is a
regular `np.ndarray`, not a `DistArray`.

When `axis` does not include the distributed dimension (axis 0), each rank
sums locally and the result is gathered. When axis 0 is included, an
`Allreduce` is performed.

---

### `DistArray.transpose` / `transpose2d`

```python
transpose(axes=None) -> DistArray
transpose2d() -> DistArray
```

Transpose the first two dimensions of the distributed array. Collective
operation. For a 2D `DistArray` with shape `(n/m, m_local)`, the result
has shape `(m_local_new, n)` where each rank holds a different row slice.

`transpose` delegates to `transpose2d` for 2D arrays. An `Alltoall` is used
internally for the optimized path; `transpose2d_ref` gathers to root for the
reference implementation.

---

### `DistArray.conj`

```python
conj() -> DistArray
```

Return a new `DistArray` with the complex conjugate of the local data.

---

## Distributed Linear Algebra

### `d_matmul`

```python
d_matmul(d_mat: DistArray, d_vec: DistArray) -> DistArray
```

Compute the distributed matrix-vector (or matrix-matrix) product
`d_mat @ d_vec`. The full `d_vec` is gathered on each rank via
`all_gather_arr`, then local `np.matmul` is performed.

| Parameter | Type | Description |
|---|---|---|
| `d_mat` | `DistArray` | Distributed matrix (first dim is distributed) |
| `d_vec` | `DistArray` | Distributed vector / matrix |
| **Returns** | `DistArray` | Result of `d_mat @ d_vec` |

---

### `d_trace`

```python
d_trace(d_mat: DistArray) -> float | np.ndarray
```

Compute the distributed trace. Each rank contributes the trace of its local
diagonal block (offset by `rank * d_vec_len`), then results are reduced via
`Allreduce`.

---

## Scatter and Gather

### `scatter_arr`

```python
scatter_arr(vec: np.ndarray, root: int = 0, comm=None) -> DistArray
```

Scatter a `np.ndarray` from `root` to all ranks, returning a `DistArray`.
The array is zero-padded if its first dimension is not evenly divisible by
the number of ranks. Only the `root` rank needs to supply valid data.

---

### `gather_arr`

```python
gather_arr(d_vec: DistArray, root: int = 0) -> np.ndarray | None
```

Gather a `DistArray` onto `root`, returning a `np.ndarray` with padded
zeros removed. Returns `None` on non-root ranks. Collective operation.

---

### `all_gather_arr`

```python
all_gather_arr(d_vec: DistArray) -> np.ndarray
```

Gather a `DistArray` onto all ranks, returning a `np.ndarray` with padded
zeros removed. Collective operation.

---

## Examples

### Basic Distributed Array Operations

```python
import numpy as np
import qlat as q
from qlat.mat_mpi import DistArray, scatter_arr, all_gather_arr, set_mpi_comm

from mpi4py import MPI

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

set_mpi_comm(MPI.COMM_WORLD)

# Scatter a vector to all ranks
vec = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
d_vec = scatter_arr(vec)
print(f"Rank {d_vec.comm.Get_rank()}: local shape = {d_vec.x.shape}")

# Arithmetic
d_sum = d_vec + d_vec
d_prod = d_vec * 2.0
full_sum = all_gather_arr(d_sum)
print(f"Gathered sum: {full_sum}")

q.end_with_mpi()
```

### Distributed Matrix-Vector Product

```python
import numpy as np
import qlat as q
from qlat.mat_mpi import DistArray, scatter_arr, all_gather_arr, set_mpi_comm

from mpi4py import MPI

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

set_mpi_comm(MPI.COMM_WORLD)

# Create a 4x4 matrix and a vector
mat = np.arange(16, dtype=np.float64).reshape(4, 4)
vec = np.ones(4, dtype=np.float64)

d_mat = scatter_arr(mat)
d_vec = scatter_arr(vec)

# Distributed matmul
d_result = d_mat @ d_vec
result = all_gather_arr(d_result)
print(f"mat @ vec = {result}")  # should match mat @ vec

q.end_with_mpi()
```

### Distributed Trace

```python
import numpy as np
import qlat as q
from qlat.mat_mpi import DistArray, scatter_arr, d_trace, set_mpi_comm

from mpi4py import MPI

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)

set_mpi_comm(MPI.COMM_WORLD)

mat = np.arange(16, dtype=np.float64).reshape(4, 4)
d_mat = scatter_arr(mat)

tr = d_trace(d_mat)
print(f"Trace = {tr}")  # should be 0 + 5 + 10 + 15 = 30

q.end_with_mpi()
```
