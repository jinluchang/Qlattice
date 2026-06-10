# `qlat.selected_points_types` — Typed Selected-Points Classes

Source: `qlat/qlat/selected_points_types.pyx.in`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Generated Types](#generated-types)
3. [Construction](#construction)
4. [Properties](#properties)
5. [Assignment and Copying](#assignment-and-copying)
6. [Mutation Operations](#mutation-operations)
7. [Arithmetic Operations](#arithmetic-operations)
8. [I/O and Serialization](#io-and-serialization)
9. [LatData Conversion](#latdata-conversion)
10. [Broadcasting](#broadcasting)
11. [NumPy Buffer Protocol](#numpy-buffer-protocol)
12. [Selected-Points Type Dictionary](#selected-points-type-dictionary)
13. [Examples](#examples)

---

## Overview

`selected_points_types` is a template module (`.pyx.in`) that generates concrete
`SelectedPoints{Type}` subclasses of `SelectedPointsBase` for every lattice data
type. Each generated class stores data at a discrete set of lattice points
identified by a `PointsSelection`.

Unlike `SelectedField` (which stores data at sites in a contiguous block
according to a `FieldSelection`), `SelectedPoints` stores data at arbitrary
lattice coordinates. This makes them suitable for correlator measurements,
point-source operators, and data exchange via shuffle plans.

All selected-points classes are registered in `selected_points_type_dict`.

## Generated Types

The template generates `SelectedPoints` variants for all 16 element types:

| Category | Classes |
|---|---|
| Complex matrices | `SelectedPointsColorMatrix`, `SelectedPointsWilsonMatrix`, `SelectedPointsNonRelWilsonMatrix`, `SelectedPointsIsospinMatrix`, `SelectedPointsSpinMatrix` |
| Complex scalars | `SelectedPointsWilsonVector`, `SelectedPointsComplexD`, `SelectedPointsComplexF` |
| Real scalars | `SelectedPointsRealD`, `SelectedPointsRealF` |
| Integer types | `SelectedPointsLong`, `SelectedPointsInt`, `SelectedPointsInt64t`, `SelectedPointsInt32t`, `SelectedPointsInt8t`, `SelectedPointsChar` |

## Construction

```python
sp = SelectedPointsRealD()                              # empty
sp = SelectedPointsRealD(n_points, multiplicity, points_dist_type)  # from sizes
sp = SelectedPointsRealD(None)                          # empty
sp = SelectedPointsRealD(psel)                          # with points selection
sp = SelectedPointsRealD(psel, multiplicity)            # with multiplicity
sp = SelectedPointsRealD(sp_src, ssp)                   # via shuffle plan
sp = SelectedPointsRealD(sf, ssp)                       # from SelectedField via shuffle
```

Parameters:

- `psel` — a `PointsSelection` specifying which points to store.
- `points_dist_type` — distribution type: `"g"` (global), `"l"` (local), `"r"` (round-robin), `"f"`, `"o"`.
- `ssp` — a `SelectedShufflePlan` for data redistribution.

## Properties

| Property | Type | Description |
|---|---|---|
| `n_points` | `int` | Number of selected points |
| `multiplicity` | `int` | Number of elements per point |
| `sizeof_m` | `int` | Size of one element in bytes |
| `points_dist_type` | `str` | Distribution type: `"g"`, `"l"`, or `"r"` |
| `geo` | `Geometry` | Geometry from `psel` |

The `points_dist_type` property is read/write.

## Assignment and Copying

The `@=` operator assigns data. It does **not** change `self.psel`.

```python
sp1 @= sp2          # SelectedPoints ← SelectedPoints (psel may differ)
sp  @= sf           # SelectedPoints ← SelectedField (gather available points)
sp  @= f            # SelectedPoints ← Field (gather at selected points)
```

| Method | Description |
|---|---|
| `copy(is_copying_data=True)` | Deep copy with same `psel` |
| `set_zero()` | Set all elements to zero |
| `swap(f1)` | Swap contents with another `SelectedPoints` of the same type |
| `swap_cast(f1)` | Swap raw bytes with a `SelectedPointsChar` |

## Mutation Operations

```python
sp.set_zero()
sp.set_rand(rng, upper=1.0, lower=0.0)        # uniform random
sp.set_rand_g(rng, center=0.0, sigma=1.0)     # Gaussian random
```

| Method | Description |
|---|---|
| `qnorm()` | Global squared norm (scalar float) |
| `qnorm_field()` | Per-point squared norm → `SelectedPointsRealD` |

## Arithmetic Operations

In-place arithmetic is supported:

```python
sp1 += sp2          # element-wise addition
sp1 -= sp2          # element-wise subtraction
sp  *= 2.0          # scalar multiplication
sp  *= (1.0+0.5j)   # complex scalar multiplication
```

## I/O and Serialization

```python
sp.save(path, is_sync_node=True)
sp.load(path, is_sync_node=True)
sp.save_str(is_sync_node=True)       # serialize to bytes (node root only)
sp.load_str(content, is_sync_node=True)  # deserialize from bytes (node 0 only)
```

## LatData Conversion

Types with double-precision, single-precision, long, or int elements support
conversion to and from `LatData`:

```python
ld = sp.to_lat_data()                # → LatData / LatDataRealF / LatDataLong / LatDataInt
sp.from_lat_data(ld)                 # ← LatData
```

The concrete `LatData` type depends on the element type:

| Element types | LatData type |
|---|---|
| `ColorMatrix`, `WilsonMatrix`, `NonRelWilsonMatrix`, `IsospinMatrix`, `SpinMatrix`, `WilsonVector`, `ComplexD`, `RealD` | `LatData` (double precision) |
| `ComplexF`, `RealF` | `LatDataRealF` (single precision) |
| `Long`, `Int64t` | `LatDataLong` |
| `Int`, `Int32t` | `LatDataInt` |

## Broadcasting

```python
sp.bcast(root=0)
```

Broadcast data from `root` to all nodes. Also broadcasts `self.psel`.
Sets `self.psel` to `None` if `psel.n_points == 0`.

## NumPy Buffer Protocol

All selected-points types support zero-copy NumPy access:

```python
arr = np.asarray(sp)    # shape: (n_points, multiplicity, *element_shape)
```

## Selected-Points Type Dictionary

```python
from qlat.field_type_dict import selected_points_type_dict

SelectedPointsRealD = selected_points_type_dict[ElemTypeRealD]
```

## Examples

```python
import qlat as q
import numpy as np

q.begin_with_mpi([[1, 1, 1, 4]])

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))

# Create a PointsSelection with a few global coordinates
total_site = q.Coordinate([4, 4, 4, 8])
psel = q.PointsSelection(total_site, [[0, 0, 0, 0], [1, 1, 1, 1]])

# Construct SelectedPoints with multiplicity 1
sp = q.SelectedPointsRealD(psel, 1)
sp.set_zero()
print(f"n_points: {sp.n_points}")

# Set uniform random values
rng = q.RngState("seed-1")
sp.set_rand(rng)

# Access as NumPy array
arr = np.asarray(sp)
print(f"Shape: {arr.shape}")

# Arithmetic
sp += sp
sp *= 0.5

# Norm
norm = sp.qnorm()
print(f"qnorm: {norm}")

# LatData round-trip
ld = sp.to_lat_data()
sp2 = q.SelectedPointsRealD(psel, 1)
sp2.from_lat_data(ld)

# Save/load
sp.save("sp_data")
sp3 = q.SelectedPointsRealD()
sp3.load("sp_data")

q.end_with_mpi()
```
