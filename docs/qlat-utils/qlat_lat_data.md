# `qlat_utils.lat_data` — Lattice Data Container

Source: `qlat-utils/qlat_utils/lat_data.pyx.in`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [LatData Types](#latdata-types)
3. [Construction](#construction)
4. [Dimension Metadata](#dimension-metadata)
5. [Data Access](#data-access)
6. [I/O Operations](#io-operations)
7. [MPI Operations](#mpi-operations)
8. [Arithmetic](#arithmetic)
9. [NumPy Buffer Protocol](#numpy-buffer-protocol)
10. [Serialization](#serialization)
11. [Examples](#examples)

---

## Overview

`lat_data` is a template module (`.pyx.in`) that generates `LatData` container
classes for multi-dimensional arrays with labeled dimensions and named indices.
Each `LatData` instance stores a ranked array where each dimension has a name
and optional string labels for its indices, making it suitable for
momentum-projected correlator data, per-site observables, and similar
structured lattice QCD datasets.

The generated classes support the NumPy buffer protocol for zero-copy array
access, MPI broadcast and global sum operations, and a human-readable
text-based file format.

---

## LatData Types

| Class | Element Type | Complex by Default |
|-------|-------------|-------------------|
| `LatData` | `RealD` (double) | Yes |
| `LatDataRealF` | `RealF` (float) | Yes |
| `LatDataInt` | `Int` (int) | No |
| `LatDataLong` | `Long` (long) | No |

All four classes share the same API. The type is selected based on the
precision and integer requirements of the data.

---

## Construction

### `mk_lat_data(info_list, is_complex=True)`

Factory function to create a `LatData` with the given dimension layout.

```python
ld = q.mk_lat_data([["t", 8], ["p", 3]])
```

| Parameter | Description |
|-----------|-------------|
| `info_list` | List of `[dim_name, dim_size, dim_indices?]` entries |
| `is_complex` | Whether the last (innermost) dimension is complex (re/im) |

### `mk_lat_data_real_f(info_list, is_complex=True)`

Same as `mk_lat_data` but creates a `LatDataRealF` (single-precision).

### `mk_lat_data_int(info_list)`

Same as `mk_lat_data` but creates a `LatDataInt` (integer, non-complex).

### `mk_lat_data_long(info_list)`

Same as `mk_lat_data` but creates a `LatDataLong` (long integer, non-complex).

### `LatData.set_info(info_list, is_complex=True)`

Set the dimension layout on an existing instance.

```python
ld = q.LatData()
ld.set_info([["t", 8], ["mom", 3, ["(0,0,0)", "(1,0,0)", "(1,1,0)"]]])
```

---

## Dimension Metadata

| Method | Description |
|--------|-------------|
| `ndim()` | Number of dimensions (excluding the implicit re/im axis) |
| `dim_name(dim)` | Name of dimension `dim` |
| `dim_size(dim)` | Size of dimension `dim` |
| `dim_names()` | List of all dimension names |
| `dim_sizes()` | List of all dimension sizes |
| `dim_indices(dim, is_filling_default=False)` | String labels for dimension `dim` |
| `dim_idx(dim, idx)` | Look up the integer index of a string label |
| `info(dim=None)` | Full info for one dimension, or all dimensions if `dim` is `None` |
| `is_complex()` | Whether the data has a complex (re/im) dimension |
| `set_dim_sizes(sizes, is_complex=True)` | Resize (data is lost) |
| `set_dim_name(dim, name, indices=None)` | Set name and optional labels for dimension `dim` |

---

## Data Access

`LatData` supports standard Python indexing via `__getitem__` and
`__setitem__`, backed by `np.asarray`:

```python
ld[0, 1] = 3.14
val = ld[0, 1]
```

### Conversion

| Method | Description |
|--------|-------------|
| `to_numpy()` | Copy data to a NumPy array |
| `from_numpy(arr, dim_names=None, is_complex=True)` | Load from a NumPy array |
| `to_list()` | Flatten to a Python list |
| `from_list(val, is_complex=True)` | Load from a flat Python list |
| `to_xarray()` | Convert to an `xarray.DataArray` with labeled coordinates |

---

## I/O Operations

| Method | Description |
|--------|-------------|
| `save(path)` | Save from node 0 |
| `load(path)` | Load from node 0 and broadcast |
| `save_node(path)` | Save on every node |
| `load_node(path)` | Load on every node (no broadcast) |
| `save_str()` | Serialize to a bytes object |
| `load_str(content)` | Deserialize from bytes |

### `load_lat_data(path)`

Factory function: load a `LatData` from file.

```python
ld = q.load_lat_data("data/prop.lat")
```

### `load_lat_data_real_f(path)`

Same but returns `LatDataRealF`.

---

## MPI Operations

| Method | Description |
|--------|-------------|
| `bcast(root=0)` | Broadcast data from `root` to all nodes |
| `glb_sum()` | Return a new `LatData` with global-summed data |
| `glb_sum_in_place()` | Global-sum in place |

---

## Arithmetic

`LatData` (the double-precision variant) supports arithmetic operations:

| Operation | Description |
|-----------|-------------|
| `ld1 + ld2` | Element-wise addition |
| `ld1 - ld2` | Element-wise subtraction |
| `ld * scalar` | Scalar multiplication (float or complex) |
| `-ld` | Negation |
| `ld += ld2` | In-place addition |
| `ld -= ld2` | In-place subtraction |
| `ld *= scalar` | In-place scalar multiplication |
| `ld.qnorm()` | Squared norm of the data vector |
| `ld.set_zero()` | Zero all elements |
| `ld.is_match(ld2)` | Check if dimensions are identical |

---

## NumPy Buffer Protocol

All `LatData` classes implement the buffer protocol, enabling zero-copy
NumPy array access:

```python
arr = np.asarray(ld)  # zero-copy view
```

The shape is `(dim0_size, dim1_size, ...)`. For complex data, the dtype is
`complex128` (or `complex64` for `LatDataRealF`) — the real/imaginary parts
are encoded in the dtype, not as an extra dimension.

---

## Serialization

`LatData` supports Python pickle via `__getstate__` / `__setstate__`. The
serialized form includes dimension metadata (names, sizes, indices) and the
data array, enabling faithful round-tripping through `pickle.dumps` /
`pickle.loads`.

---

## Examples

### Create and populate a LatData

```python
import qlat_utils as q

ld = q.mk_lat_data([["t", 8], ["mom", 3]], is_complex=True)
ld.set_dim_name(1, "mom", ["(0,0,0)", "(1,0,0)", "(1,1,0)"])
ld[0, 0] = 1.0 + 0.5j
if q.get_id_node() == 0:
    print("shape:", ld.dim_sizes())
    print("value:", ld[0, 0])
```

### Query dimension metadata

```python
import qlat_utils as q

ld = q.mk_lat_data([["t", 4], ["mom", 3, ["p0", "p1", "p2"]]], is_complex=True)
if q.get_id_node() == 0:
    print("ndim:", ld.ndim())               # 2
    print("dim_names:", ld.dim_names())     # ['t', 'mom']
    print("dim_sizes:", ld.dim_sizes())     # [4, 3]
    print("is_complex:", ld.is_complex())   # True
    print("dim_indices(1):", ld.dim_indices(1))  # ['p0', 'p1', 'p2']
    print("dim_idx(1, 'p1'):", ld.dim_idx(1, "p1"))  # 1
```

### Arithmetic operations

```python
import qlat_utils as q

ld_a = q.mk_lat_data([["t", 4]])
ld_b = q.mk_lat_data([["t", 4]])
ld_a[0] = 1.0
ld_b[0] = 2.0

ld_sum = ld_a + ld_b         # element-wise addition
ld_diff = ld_a - ld_b        # element-wise subtraction
ld_scaled = ld_a * 3.0       # scalar multiplication
ld_neg = -ld_a               # negation

ld_a += ld_b                 # in-place addition
ld_a *= 0.5                  # in-place scalar multiplication
if q.get_id_node() == 0:
    print("qnorm:", ld_a.qnorm())  # squared norm
```

### NumPy buffer protocol (zero-copy)

```python
import qlat_utils as q
import numpy as np

ld = q.mk_lat_data([["t", 4], ["x", 3]])
ld[0, 0] = 5.0
arr = np.asarray(ld)    # zero-copy view
arr[1, 1] = 7.0         # writes through to LatData
if q.get_id_node() == 0:
    print("ld[1,1]:", ld[1, 1])  # 7.0
```

### Pickle round-trip

```python
import qlat_utils as q
import pickle

ld = q.mk_lat_data([["t", 4], ["mom", 2, ["p0", "p1"]]])
ld[0, 0] = 1.0 + 0.5j
data = pickle.dumps(ld)
ld_restored = pickle.loads(data)
if q.get_id_node() == 0:
    print("dim_names:", ld_restored.dim_names())
    print("dim_indices(1):", ld_restored.dim_indices(1))
```

### save_str / load_str round-trip

```python
import qlat_utils as q

ld = q.mk_lat_data([["t", 4], ["p", 2]])
ld[0, 0] = 1.5
content = ld.save_str()     # serialize to bytes
ld2 = q.LatData()
ld2.load_str(content)       # deserialize
```

### Load from file and convert to NumPy

```python
import qlat_utils as q
import numpy as np

ld = q.load_lat_data("data/pion_corr.lat")
arr = np.asarray(ld)
if q.get_id_node() == 0:
    print("NumPy shape:", arr.shape)
```

### MPI global sum

```python
import qlat_utils as q

size_node_list = [[1, 1, 1, 1]]
q.begin_with_mpi(size_node_list)
ld = q.mk_lat_data([["t", 8]])
ld.set_zero()
# ... accumulate local data ...
ld.glb_sum_in_place()  # sum across all nodes
q.end_with_mpi()
```
