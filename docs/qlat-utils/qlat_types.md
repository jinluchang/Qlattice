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

All subclasses provide the following class-level attributes:

| Attribute | Type | Description |
|---|---|---|
| `name` | `str` | Human-readable type name |
| `sizeof_m` | `int` | Total byte size of the C++ type |

All subclasses also provide the following Python-accessible static methods (wrappers around the internal `cdef` methods):

| Method | Return Type | Description |
|---|---|---|
| `py_format()` | `str` | Buffer protocol format string (e.g., `'Zd'` for complex double) |
| `py_itemsize()` | `int` | Byte size of one scalar element |
| `py_ndim()` | `int` | Number of dimensions |
| `py_shape()` | `tuple` | Dimension sizes |
| `py_size()` | `int` | Total byte size of the type |

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

> **Note:** The tables below document the full type metadata. `name` and `sizeof_m` are class attributes. `py_format()`, `py_itemsize()`, `py_ndim()`, `py_shape()`, and `py_size()` are Python-accessible static methods (wrappers around internal `cdef` methods).

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

# ColorMatrix: 3x3 complex double, total size 144 bytes
print(ElemTypeColorMatrix.name)       # "ColorMatrix"
print(ElemTypeColorMatrix.sizeof_m)   # 144

# WilsonMatrix: 12x12 complex double
print(ElemTypeWilsonMatrix.name)      # "WilsonMatrix"
print(ElemTypeWilsonMatrix.sizeof_m)  # 2304

# Scalar type
print(ElemTypeRealD.name)            # "RealD"
print(ElemTypeRealD.sizeof_m)        # 8
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
# ColorMatrix is 3x3 complex double (sizeof_m = 144 bytes)
shape = (3, 3)
dtype = np.complex128                 # matches 'Zd' format
mat = np.zeros(shape, dtype=dtype)
print(mat.nbytes)                     # 144, matches ElemTypeColorMatrix.sizeof_m
```

### Querying Type Properties via Python-Accessible Methods

All ElemType subclasses expose `py_format()`, `py_itemsize()`, `py_ndim()`,
`py_shape()`, and `py_size()` as static methods callable from Python:

```python
from qlat_utils.types import ElemTypeColorMatrix, ElemTypeWilsonMatrix
from qlat_utils.types import ElemTypeComplexD, ElemTypeRealD, ElemTypeLong

# ColorMatrix: 3x3 complex double
assert ElemTypeColorMatrix.py_ndim() == 2
assert ElemTypeColorMatrix.py_shape() == (3, 3)
assert ElemTypeColorMatrix.py_itemsize() == 16
assert ElemTypeColorMatrix.py_size() == 144
assert ElemTypeColorMatrix.py_format() == 'Zd'

# WilsonMatrix: 12x12 complex double
assert ElemTypeWilsonMatrix.py_ndim() == 2
assert ElemTypeWilsonMatrix.py_shape() == (12, 12)
assert ElemTypeWilsonMatrix.py_itemsize() == 16
assert ElemTypeWilsonMatrix.py_size() == 2304

# Scalar types have ndim == 0 and empty shape
assert ElemTypeComplexD.py_ndim() == 0
assert ElemTypeComplexD.py_shape() == ()
assert ElemTypeComplexD.py_itemsize() == 16
assert ElemTypeComplexD.py_size() == 16

assert ElemTypeRealD.py_ndim() == 0
assert ElemTypeRealD.py_format() == 'd'
assert ElemTypeRealD.py_itemsize() == 8

assert ElemTypeLong.py_ndim() == 0
assert ElemTypeLong.py_format() == 'q'
assert ElemTypeLong.py_itemsize() == 8
```

### Building a NumPy Array from Type Metadata

Use the Python-accessible methods to programmatically construct arrays whose
layout matches a given ElemType, without hardcoding shape or dtype:

```python
import numpy as np
from qlat_utils.types import ElemTypeSpinMatrix, ElemTypeWilsonVector

fmt_to_dtype = {
    'Zd': np.complex128,
    'Zf': np.complex64,
    'd': np.float64,
    'f': np.float32,
    'q': np.int64,
    'i': np.int32,
    'b': np.int8,
}

def make_array(elem_type, count):
    dtype = fmt_to_dtype[elem_type.py_format()]
    shape = (count,) + elem_type.py_shape()
    return np.zeros(shape, dtype=dtype)

# 10 SpinMatrix elements: shape (10, 4, 4), dtype complex128
arr = make_array(ElemTypeSpinMatrix, 10)
assert arr.shape == (10, 4, 4)
assert arr.dtype == np.complex128
assert arr.nbytes == 10 * ElemTypeSpinMatrix.py_size()

# 5 WilsonVector elements: shape (5, 12), dtype complex128
vec = make_array(ElemTypeWilsonVector, 5)
assert vec.shape == (5, 12)
assert vec.nbytes == 5 * ElemTypeWilsonVector.py_size()
```

### Iterating Over All ElemType Subclasses

```python
import numpy as np
from qlat_utils.types import (
    ElemTypeColorMatrix, ElemTypeWilsonMatrix, ElemTypeNonRelWilsonMatrix,
    ElemTypeIsospinMatrix, ElemTypeSpinMatrix, ElemTypeWilsonVector,
    ElemTypeComplexD, ElemTypeComplexF, ElemTypeRealD, ElemTypeRealF,
    ElemTypeLong, ElemTypeInt, ElemTypeChar, ElemTypeInt64t,
    ElemTypeInt32t, ElemTypeInt8t,
)

all_types = [
    ElemTypeColorMatrix, ElemTypeWilsonMatrix, ElemTypeNonRelWilsonMatrix,
    ElemTypeIsospinMatrix, ElemTypeSpinMatrix, ElemTypeWilsonVector,
    ElemTypeComplexD, ElemTypeComplexF, ElemTypeRealD, ElemTypeRealF,
    ElemTypeLong, ElemTypeInt, ElemTypeChar, ElemTypeInt64t,
    ElemTypeInt32t, ElemTypeInt8t,
]

for t in all_types:
    assert t.py_size() == t.sizeof_m, f"{t.name}: py_size()={t.py_size()} != sizeof_m={t.sizeof_m}"
    if t.py_ndim() > 0:
        expected = t.py_itemsize() * int(np.prod(t.py_shape()))
    else:
        expected = t.py_itemsize()
    assert t.py_size() == expected, f"{t.name}: py_size()={t.py_size()} != expected={expected}"
```

### Verifying Consistency Between Properties

```python
import numpy as np
from qlat_utils.types import (
    ElemTypeColorMatrix, ElemTypeWilsonMatrix, ElemTypeNonRelWilsonMatrix,
    ElemTypeIsospinMatrix, ElemTypeSpinMatrix, ElemTypeWilsonVector,
    ElemTypeComplexD, ElemTypeComplexF, ElemTypeRealD, ElemTypeRealF,
    ElemTypeLong, ElemTypeInt, ElemTypeChar, ElemTypeInt64t,
    ElemTypeInt32t, ElemTypeInt8t,
)

all_types = [
    ElemTypeColorMatrix, ElemTypeWilsonMatrix, ElemTypeNonRelWilsonMatrix,
    ElemTypeIsospinMatrix, ElemTypeSpinMatrix, ElemTypeWilsonVector,
    ElemTypeComplexD, ElemTypeComplexF, ElemTypeRealD, ElemTypeRealF,
    ElemTypeLong, ElemTypeInt, ElemTypeChar, ElemTypeInt64t,
    ElemTypeInt32t, ElemTypeInt8t,
]

for t in all_types:
    assert t.py_size() == t.sizeof_m, f"{t.name}: py_size()={t.py_size()} != sizeof_m={t.sizeof_m}"
    if t.py_ndim() > 0:
        expected = t.py_itemsize() * int(np.prod(t.py_shape()))
    else:
        expected = t.py_itemsize()
    assert t.py_size() == expected, f"{t.name}: py_size()={t.py_size()} != expected={expected}"
    assert isinstance(t.name, str) and len(t.name) > 0
    assert isinstance(t.py_ndim(), int) and t.py_ndim() >= 0
    assert isinstance(t.py_shape(), tuple)
    assert isinstance(t.py_itemsize(), int) and t.py_itemsize() > 0
    assert isinstance(t.py_size(), int) and t.py_size() > 0
    assert isinstance(t.py_format(), str) and len(t.py_format()) > 0
```
