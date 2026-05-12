# `qlat.field_types` — Typed Lattice Field Classes

Source: `qlat/qlat/field_types.pyx.in`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Field Types](#field-types)
   - [Complex Matrix Types](#complex-matrix-types)
   - [Complex Scalar Types](#complex-scalar-types)
   - [Real Scalar Types](#real-scalar-types)
   - [Integer Types](#integer-types)
3. [Construction](#construction)
4. [Properties](#properties)
5. [Mutation Operations](#mutation-operations)
6. [Global Reductions](#global-reductions)
7. [Element Access by Coordinate](#element-access-by-coordinate)
8. [Field Shifting and Reflection](#field-shifting-and-reflection)
9. [I/O Operations](#io-operations)
10. [Assignment and Copying](#assignment-and-copying)
11. [NumPy Buffer Protocol](#numpy-buffer-protocol)
12. [Field Type Dictionary](#field-type-dictionary)
13. [Examples](#examples)

---

## Overview

`field_types` is a template module (`.pyx.in`) that generates concrete
`Field{Type}` subclasses of `FieldBase` for every lattice data type. Each
generated class (e.g., `FieldColorMatrix`, `FieldRealD`, `FieldChar`) shares a
uniform API and supports the NumPy buffer protocol for zero-copy array access.

All field classes are registered in the `field_type_dict` lookup:
`field_type_dict[ElemTypeRealD]` returns `FieldRealD`, enabling generic
construction via `Field(ctype, geo, multiplicity)`.

## Field Types

### Complex Matrix Types

| Class | Element | Element Shape |
|---|---|---|
| `FieldColorMatrix` | 3×3 complex (SU(3)) | `(2, 3)` |
| `FieldWilsonMatrix` | 4×3 complex (spin-color) | `(2, 12)` |
| `FieldNonRelWilsonMatrix` | 2×3 complex | `(2, 6)` |
| `FieldIsospinMatrix` | 2×2 complex | `(2, 2)` |
| `FieldSpinMatrix` | 4×4 complex | `(2, 4)` |

### Complex Scalar Types

| Class | Element | Element Shape |
|---|---|---|
| `FieldWilsonVector` | Complex vector, length 12 | `(1, 12)` |
| `FieldComplexD` | Double-precision complex scalar | scalar |
| `FieldComplexF` | Single-precision complex scalar | scalar |

### Real Scalar Types

| Class | Element | Element Shape |
|---|---|---|
| `FieldRealD` | Double-precision real scalar | scalar |
| `FieldRealF` | Single-precision real scalar | scalar |

### Integer Types

| Class | Element |
|---|---|
| `FieldLong` | 64-bit signed integer |
| `FieldInt` | 32-bit signed integer |
| `FieldInt64t` | 64-bit signed integer |
| `FieldInt32t` | 32-bit signed integer |
| `FieldInt8t` | 8-bit signed integer |
| `FieldChar` | 8-bit signed integer |

## Construction

```python
f = FieldRealD()                        # empty (uninitialized)
f = FieldRealD(geo)                     # single multiplicity
f = FieldRealD(geo, multiplicity=3)     # multiple multiplicities
```

`geo` is a `Geometry` object defining the lattice on which the field lives.
`multiplicity` specifies how many elements per site (default 1).

## Properties

| Property | Type | Description |
|---|---|---|
| `geo` | `Geometry` | The lattice geometry of the field |
| `total_site` | `Coordinate` | Total lattice dimensions |
| `n_sites` | `int` | Number of local sites stored on this node |
| `multiplicity` | `int` | Number of elements per site |
| `sizeof_m` | `int` | Size of one element in bytes |

## Mutation Operations

```python
f.set_zero()                            # set all elements to zero
f.set_unit(coef=1.0+0j)                 # set to unit matrix/scalar (times coef)
f.set_rand(rng, upper=1.0, lower=0.0)   # uniform random in [lower, upper)
f.set_rand_g(rng, center=0.0, sigma=1.0)  # Gaussian random
```

`rng` is an `RngState` object for reproducible random number generation.

## Global Reductions

These methods perform MPI reductions across all nodes:

```python
sp = f.glb_sum()             # global sum → SelectedPoints (n_points=1)
sp = f.glb_sum_tslice(t_dir=3)  # sum per timeslice → SelectedPoints

# Available only for RealD, RealF, Long, Int, Int64t, Int32t:
sp = f.glb_max()             # global maximum → SelectedPoints (n_points=1)
sp = f.glb_min()             # global minimum → SelectedPoints (n_points=1)
```

Each returns a `SelectedPoints{Type}` with the result.

```python
norm = f.qnorm()             # global squared norm (scalar float)
f_n = f.qnorm_field()        # per-site squared norm → FieldRealD
```

## Element Access by Coordinate

Read/write individual elements or collections of elements specified by
coordinate arrays:

```python
sp = f.get_elems_xg(xg_arr)       # get elements at coordinates → SelectedPoints
f.set_elems_xg(xg_arr, val)       # set elements at coordinates from value

sp = f.get_elem_xg(xg_arr, m)     # get single multiplicity m at coordinates
f.set_elem_xg(xg_arr, m, val)     # set single multiplicity m at coordinates
```

`xg_arr` can be an `xg` coordinate, `xg_list`, or `xg_arr`.

## Field Shifting and Reflection

```python
f1 = f.shift(shift, is_reflect=False)  # return shifted copy (self is unchanged)
```

If `shift` is `None`, returns a copy. Otherwise, the field is shifted
cyclically by the given coordinate. If `is_reflect=True`, coordinates are
negated after shifting.

## I/O Operations

```python
f.read_direct(filename, new_size_node=None)
f.write_direct(filename, new_size_node=None)

f.read_sfr_dynamic(sfr, filename)           # read from ShuffledFieldsReader
f.write_sfw_dynamic(sfw, filename, skip_if_exist=False)  # write via ShuffledFieldsWriter
```

## Assignment and Copying

The `@=` operator assigns field data. It accepts `Field`, `SelectedField`, or
`SelectedPoints` of the same type:

```python
f1 @= f2                         # copy field data (geo unchanged if already set)
f1 @= selected_field              # scatter SelectedField into Field
f1 @= selected_points             # scatter SelectedPoints into Field
```

Other assignment methods:

```python
f2 = f.copy(is_copying_data=True)  # deep copy (or shallow if False)
f.swap(f1)                         # swap data with another Field{Type}
f.swap_cast(f1)                    # swap with FieldChar (byte-cast)
f.swap_sp_cast(sp, geo)            # swap with SelectedPointsChar (byte-cast)
```

## NumPy Buffer Protocol

All field types implement the Python buffer protocol, enabling zero-copy NumPy
access:

```python
arr = np.asarray(f)               # shape: (n_sites, multiplicity, *element_shape)
arr = np.asarray(f, dtype=np.float64)       # real components only
arr = np.asarray(f, dtype=np.complex128)    # complex components only
```

For scalar types (RealD, ComplexD, Int, etc.), the element shape is empty.
For matrix types, element dimensions appear after the multiplicity axis.

## Field Type Dictionary

The module registers each generated class into `field_type_dict`:

```python
from qlat.field_type_dict import field_type_dict

FieldRealD = field_type_dict[ElemTypeRealD]     # equivalent to FieldRealD
f = Field(ElemTypeRealD, geo, multiplicity)     # generic construction via factory
```

This enables generic code that constructs fields by element type.

## Examples

```python
import qlat as q
import numpy as np

q.begin_with_mpi([1, 1, 1, 4])

# Create a geometry and field
geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))
f = q.FieldRealD(geo)

# Set to Gaussian random
rng = q.RngState()
rng.reset(42)
f.set_rand_g(rng, center=0.0, sigma=1.0)

# Access as NumPy array (zero-copy)
arr = np.asarray(f)
print(f"Shape: {arr.shape}")            # (local_sites,)

# Set zero and set a site to 1
f.set_zero()
f.set_elem_xg([0, 0, 0, 0], 0, 1.0)
sp = f.get_elem_xg([0, 0, 0, 0], 0)
print(f"Value: {np.asarray(sp)}")

# Global sum
sp_sum = f.glb_sum()
total = np.asarray(sp_sum).ravel()[0]
print(f"Total sum: {total}")

# Copy and shift
f_shift = f.shift(q.Coordinate([1, 0, 0, 0]))

# I/O
f.write_direct("field.dat")
f2 = q.FieldRealD()
f2.read_direct("field.dat")

q.end_with_mpi()
```
