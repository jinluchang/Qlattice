# `qlat_utils.types` — Element Type Descriptors and Buffer Protocol Support

Source: `qlat-utils/qlat_utils/types.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Buffer Class](#buffer-class)
3. [ElemType Hierarchy](#elemtype-hierarchy)
4. [Matrix Types](#matrix-types)
5. [Scalar Types](#scalar-types)
6. [Type Properties Reference](#type-properties-reference)
7. [Examples](#examples)

---

## Overview

`types` defines the Python-level type descriptors for the fundamental data types used in lattice QCD calculations. It provides:

- **`Buffer`** — a Cython-level implementation of the Python buffer protocol (`__getbuffer__` / `__releasebuffer__`) that allows C++ memory to be exposed as NumPy-compatible buffers without copying.
- **`ElemType` class hierarchy** — a family of descriptor classes, each encoding the shape, itemsize, format string, and byte-size of a specific C++ lattice type (e.g., `ColorMatrix`, `WilsonMatrix`, `ComplexD`, `RealD`, etc.).

These types serve as metadata bridges between C++ field storage and Python's buffer/NumPy ecosystem.

Python access:

```python
import qlat_utils as q
```

---

## Buffer Class

```cython
cdef class Buffer:
    def __cinit__(self, object obj=None, int ndim=1, int itemsize=1,
                  char* fmt=NULL, char* buf=NULL)
```

A Cython extension type that implements the Python buffer protocol. Instances are not typically created directly by users; they are used internally by field types to expose their data as buffer-compatible objects.

### Constructor Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `obj` | `object` | `None` | Parent object that owns the memory. Its `release_buffer()` method is called when the buffer is released. |
| `ndim` | `int` | `1` | Number of dimensions |
| `itemsize` | `int` | `1` | Size in bytes of a single element |
| `fmt` | `char*` | `NULL` | Format string for the buffer protocol (e.g. `'Zd'` for complex double) |
| `buf` | `char*` | `NULL` | Pointer to the raw memory |

### Internal Methods

| Method | Description |
|---|---|
| `set_buffer(buffer, flags)` | Populate a `Py_buffer` struct for the buffer protocol |
| `set_dim_size(dim, size)` | Set the size of dimension `dim` |
| `update_strides_from_shape()` | Compute C-contiguous strides from the current shape |
| `get_len()` | Return total byte length (`itemsize * product of shape`) |

### Buffer Protocol

`Buffer` implements `__getbuffer__` and `__releasebuffer__`, making it compatible with `memoryview`, `np.asarray()`, and any Python API that accepts the buffer protocol. When the buffer is released, the parent object's `release_buffer()` method is called.

---

## ElemType Hierarchy

`ElemType` is the base class for all element type descriptors. Each subclass encodes the metadata for a specific C++ lattice data type.

```python
class ElemType:
    name = ""
```

All subclasses provide the following class-level attributes and static methods:

| Attribute/Method | Type | Description |
|---|---|---|
| `name` | `str` | Human-readable type name |
| `sizeof_m` | `int` | Total byte size of the C++ type |
| `format()` | `str` | Buffer protocol format string |
| `itemsize()` | `int` | Byte size of one scalar element |
| `ndim()` | `int` | Number of dimensions when viewed as an array of scalars |
| `shape()` | `tuple` | Shape when viewed as an array of scalars |
| `size()` | `int` | Total byte size (same as `sizeof_m`) |

---

## Matrix Types

These types represent 2D matrix structures used in lattice QCD calculations. All use `ComplexD` (complex double, `'Zd'`) as their scalar element type.

### `ElemTypeColorMatrix`

3×3 complex color matrix. Used for SU(3) gauge links.

| Property | Value |
|---|---|
| `name` | `"ColorMatrix"` |
| `ndim` | `2` |
| `shape` | `(3, 3)` |
| `itemsize` | `sizeof(ComplexD)` = 16 bytes |
| `size` | 144 bytes |

### `ElemTypeWilsonMatrix`

12×12 complex Wilson (Dirac-spinor × color) matrix.

| Property | Value |
|---|---|
| `name` | `"WilsonMatrix"` |
| `ndim` | `2` |
| `shape` | `(12, 12)` |
| `itemsize` | `sizeof(ComplexD)` = 16 bytes |
| `size` | 2304 bytes |

### `ElemTypeNonRelWilsonMatrix`

6×6 complex non-relativistic Wilson matrix.

| Property | Value |
|---|---|
| `name` | `"NonRelWilsonMatrix"` |
| `ndim` | `2` |
| `shape` | `(6, 6)` |
| `itemsize` | `sizeof(ComplexD)` = 16 bytes |
| `size` | 576 bytes |

### `ElemTypeIsospinMatrix`

2×2 complex isospin matrix.

| Property | Value |
|---|---|
| `name` | `"IsospinMatrix"` |
| `ndim` | `2` |
| `shape` | `(2, 2)` |
| `itemsize` | `sizeof(ComplexD)` = 16 bytes |
| `size` | 64 bytes |

### `ElemTypeSpinMatrix`

4×4 complex Dirac spin matrix.

| Property | Value |
|---|---|
| `name` | `"SpinMatrix"` |
| `ndim` | `2` |
| `shape` | `(4, 4)` |
| `itemsize` | `sizeof(ComplexD)` = 16 bytes |
| `size` | 256 bytes |

### `ElemTypeWilsonVector`

12-component complex Wilson vector (Dirac spin × color).

| Property | Value |
|---|---|
| `name` | `"WilsonVector"` |
| `ndim` | `1` |
| `shape` | `(12,)` |
| `itemsize` | `sizeof(ComplexD)` = 16 bytes |
| `size` | 192 bytes |

---

## Scalar Types

These types represent scalar (0-dimensional) primitive types.

### `ElemTypeComplexD`

Complex double precision (`std::complex<double>`).

| Property | Value |
|---|---|
| `name` | `"ComplexD"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'Zd'` |
| `itemsize` | 16 bytes |

### `ElemTypeComplexF`

Complex single precision (`std::complex<float>`).

| Property | Value |
|---|---|
| `name` | `"ComplexF"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'Zf'` |
| `itemsize` | 8 bytes |

### `ElemTypeRealD`

Double precision real (`double`).

| Property | Value |
|---|---|
| `name` | `"RealD"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'d'` |
| `itemsize` | 8 bytes |

### `ElemTypeRealF`

Single precision real (`float`).

| Property | Value |
|---|---|
| `name` | `"RealF"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'f'` |
| `itemsize` | 4 bytes |

### `ElemTypeLong`

64-bit signed integer (`Long`).

| Property | Value |
|---|---|
| `name` | `"Long"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'q'` |
| `itemsize` | 8 bytes |

### `ElemTypeInt`

32-bit signed integer (`int`).

| Property | Value |
|---|---|
| `name` | `"Int"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'i'` |
| `itemsize` | 4 bytes |

### `ElemTypeChar`

Signed character (`char`).

| Property | Value |
|---|---|
| `name` | `"Char"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'b'` |
| `itemsize` | 1 byte |

### `ElemTypeInt64t`

64-bit signed integer (`int64_t`).

| Property | Value |
|---|---|
| `name` | `"Int64t"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'q'` |
| `itemsize` | 8 bytes |

### `ElemTypeInt32t`

32-bit signed integer (`int32_t`).

| Property | Value |
|---|---|
| `name` | `"Int32t"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'i'` |
| `itemsize` | 4 bytes |

### `ElemTypeInt8t`

8-bit signed integer (`int8_t`).

| Property | Value |
|---|---|
| `name` | `"Int8t"` |
| `ndim` | `0` |
| `shape` | `()` |
| `format` | `'b'` |
| `itemsize` | 1 byte |

---

## Type Properties Reference

### Matrix Types Summary

| Type | Name | Shape | Item Format | Element Size | Total Size |
|---|---|---|---|---|---|
| `ElemTypeColorMatrix` | `ColorMatrix` | `(3, 3)` | `Zd` | 16 B | 144 B |
| `ElemTypeWilsonMatrix` | `WilsonMatrix` | `(12, 12)` | `Zd` | 16 B | 2304 B |
| `ElemTypeNonRelWilsonMatrix` | `NonRelWilsonMatrix` | `(6, 6)` | `Zd` | 16 B | 576 B |
| `ElemTypeIsospinMatrix` | `IsospinMatrix` | `(2, 2)` | `Zd` | 16 B | 64 B |
| `ElemTypeSpinMatrix` | `SpinMatrix` | `(4, 4)` | `Zd` | 16 B | 256 B |
| `ElemTypeWilsonVector` | `WilsonVector` | `(12,)` | `Zd` | 16 B | 192 B |

### Scalar Types Summary

| Type | Name | Format | Size |
|---|---|---|---|
| `ElemTypeComplexD` | `ComplexD` | `Zd` | 16 B |
| `ElemTypeComplexF` | `ComplexF` | `Zf` | 8 B |
| `ElemTypeRealD` | `RealD` | `d` | 8 B |
| `ElemTypeRealF` | `RealF` | `f` | 4 B |
| `ElemTypeLong` | `Long` | `q` | 8 B |
| `ElemTypeInt` | `Int` | `i` | 4 B |
| `ElemTypeChar` | `Char` | `b` | 1 B |
| `ElemTypeInt64t` | `Int64t` | `q` | 8 B |
| `ElemTypeInt32t` | `Int32t` | `i` | 4 B |
| `ElemTypeInt8t` | `Int8t` | `b` | 1 B |

---

## Examples

### Querying Type Metadata

```python
import qlat_utils as q
from qlat_utils.types import ElemTypeColorMatrix, ElemTypeWilsonMatrix, ElemTypeRealD

# ColorMatrix: 3x3 complex double
print(ElemTypeColorMatrix.name)       # "ColorMatrix"
print(ElemTypeColorMatrix.ndim())     # 2
print(ElemTypeColorMatrix.shape())    # (3, 3)
print(ElemTypeColorMatrix.itemsize()) # 16  (sizeof complex double)
print(ElemTypeColorMatrix.size())     # 144 (16 * 9)

# WilsonMatrix: 12x12 complex double
print(ElemTypeWilsonMatrix.shape())   # (12, 12)
print(ElemTypeWilsonMatrix.size())    # 2304

# Scalar type
print(ElemTypeRealD.name)            # "RealD"
print(ElemTypeRealD.ndim())          # 0
print(ElemTypeRealD.format())        # 'd'
```

### Buffer Protocol Usage

The `Buffer` class is used internally by field types. A typical pattern:

```python
import numpy as np
import qlat_utils as q
from qlat_utils.types import Buffer, ElemTypeComplexD

# Buffer is typically created internally by field objects.
# When a field exposes its data via the buffer protocol,
# NumPy can access it without copying:
#   arr = np.asarray(field_buffer)
# The parent object controls memory lifetime via release_buffer().
```

### Using ElemType for Array Construction

```python
import numpy as np
from qlat_utils.types import ElemTypeColorMatrix

# Construct a zero-filled array matching the ColorMatrix layout
shape = ElemTypeColorMatrix.shape()   # (3, 3)
dtype = np.complex128                 # matches 'Zd' format
mat = np.zeros(shape, dtype=dtype)
print(mat.nbytes)                     # 144, matches ElemTypeColorMatrix.size()
```
