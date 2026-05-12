# `qlat.field_base` — Base Classes and Factories for Lattice Fields

Source: `qlat/qlat/field_base.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Factory Functions](#factory-functions)
   - [`Field`](#field)
   - [`SelectedField`](#selectedfield)
   - [`SelectedPoints`](#selectedpoints)
3. [`FieldBase` Class](#fieldbase-class)
   - [Type Casting](#type-casting)
   - [Arithmetic Operators](#arithmetic-operators)
   - [Element Access](#element-access)
   - [Serialization — Save and Load](#serialization--save-and-load)
   - [Endianness and Precision Conversion](#endianness-and-precision-conversion)
   - [Checksums](#checksums)
   - [Pickle Support](#pickle-support)
4. [`SelectedFieldBase` Class](#selectedfieldbase-class)
5. [`SelectedPointsBase` Class](#selectedpointsbase-class)
6. [Field Manipulation Utilities](#field-manipulation-utilities)
   - [`split_fields`](#split_fields)
   - [`merge_fields`](#merge_fields)
   - [`merge_fields_ms`](#merge_fields_ms)
   - [`mk_merged_fields_ms`](#mk_merged_fields_ms)
7. [Examples](#examples)

---

## Overview

`field_base` defines the three core abstract base classes and their
corresponding factory functions that form the data layer of qlat:

| Class | Description | Factory |
|---|---|---|
| `FieldBase` | A lattice field defined on every site of a `Geometry` | `Field(ctype, geo, multiplicity)` |
| `SelectedFieldBase` | A field defined only on sites selected by a `FieldSelection` | `SelectedField(ctype, fsel, multiplicity)` |
| `SelectedPointsBase` | Data at a discrete set of points selected by a `PointsSelection` | `SelectedPoints(ctype, psel, multiplicity)` |

Concrete subclasses (e.g., `FieldRealD`, `FieldComplexD`, `SelectedFieldRealF`,
`SelectedPointsComplexD`, …) are defined in `field_types.pyx`,
`selected_field_types.pyx`, and `selected_points_types.pyx`. The factory
functions here provide a uniform way to construct them by `ctype`.

Each field stores an array of elements indexed by site (and a multiplicity
index within each site). Elements can be real/complex scalars or small fixed-size
matrices (e.g., `ColorMatrix`, `WilsonMatrix`). The `ctype` determines the
element type.

All field classes support NumPy array interface (`np.asarray(field)`) for
zero-copy interoperability with NumPy operations.

---

## Factory Functions

### `Field`

```python
Field(ctype, geo=None, multiplicity=0) -> FieldBase subclass
```

Create a `Field` of the given element type.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ctype` | `ElemType` | — | Element type (e.g., `ElemTypeRealD`, `ElemTypeComplexD`) |
| `geo` | `Geometry` or `None` | `None` | Lattice geometry; if `None`, field is uninitialized |
| `multiplicity` | `int` | `0` | Number of elements per site |

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
f = q.Field(q.ElemTypeRealD, geo, 4)  # 4 real doubles per site

q.end_with_mpi()
```

### `SelectedField`

```python
SelectedField(ctype, fsel, multiplicity=0) -> SelectedFieldBase subclass
```

Create a `SelectedField` defined only on the sites selected by `fsel`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ctype` | `ElemType` | — | Element type |
| `fsel` | `FieldSelection` | — | Site selection |
| `multiplicity` | `int` | `0` | Number of elements per selected site |

### `SelectedPoints`

```python
SelectedPoints(ctype, psel, multiplicity=0) -> SelectedPointsBase subclass
```

Create a `SelectedPoints` holding data at the points specified by `psel`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ctype` | `ElemType` | — | Element type |
| `psel` | `PointsSelection` | — | Point selection |
| `multiplicity` | `int` | `0` | Number of elements per point |

---

## `FieldBase` Class

`FieldBase` is the abstract base class for all concrete `Field` types. It
provides common operations shared by all element types.

### Type Casting

#### `cast_from(other: FieldBase)`

Cast the raw bytes of `other` (which may have a different element type) into
`self`. The total byte size per site must match:
`other.multiplicity * other.sizeof_m == self.multiplicity * self.sizeof_m`.

`self` is re-initialized with `other.geo` and the appropriate multiplicity.

#### `get_data_sig(rng: RngState) -> float | complex`

Compute a deterministic signature of the field data by dotting with a random
±1 vector drawn from `rng`. Casts internally to `ComplexD` or `RealD` as
needed. The result is globally summed across all MPI ranks.

Useful for checksumming field data to detect corruption or divergence.

#### `as_field(ctype=ElemTypeComplexD) -> FieldBase`

Return a new `Field` of the specified `ctype` with the same content. Performs
type conversion if needed.

#### `from_field(f)`

Assign content from `f` into `self`, converting types if needed.

### Arithmetic Operators

`FieldBase` supports in-place arithmetic with other `FieldBase`,
`SelectedFieldBase`, and `SelectedPointsBase` objects.

#### `__iadd__` (`+=`)

```python
field += other_field
field += selected_field    # accumulates only at selected sites
field += selected_points   # accumulates at the corresponding sites
```

#### `__isub__` (`-=`)

```python
field -= other_field
field -= selected_field
field -= selected_points
```

#### `__imul__` (`*=`)

```python
field *= 2.0                        # scalar multiplication
field *= complex(1.0, 0.5)          # complex scalar multiplication
field *= other_field                # element-wise (other must be ComplexD or RealD)
```

### Element Access

Elements are accessed via NumPy array view (`np.asarray(self)`).

#### `__getitem__(idx)` / `__setitem__(idx, val)`

Index into the field as a NumPy array. Typical shapes:
- `(n_sites, multiplicity, ...)` where `...` depends on the element type.

```python
arr = f[:]            # full array view
val = f[0, 0]         # site 0, multiplicity 0
f[0, 0] = np.array([1.0, 2.0, 3.0])
```

#### `get_elems(idx)` / `set_elems(idx, val)`

Get/set all multiplicity elements at site index `idx`.

#### `get_elem(idx, m=0)` / `set_elem(idx, m, val)`

Get/set a single element at site `idx`, multiplicity `m`.

#### `set_m(f1, m, m1)`

Set component `m` of `self` from component `m1` of `f1`.

### Serialization — Save and Load

Fields support direct (raw binary) and endianness-aware save/load, both to
local files and to shuffled field archives.

All `save_*` and `load_*` methods accept two calling conventions:
- **Local file:** `f.save_64("path/to/file")`
- **Shuffled archive:** `f.save_64(sfw, "field_name")` / `f.load_64(sfr, "field_name")`

#### `save_direct(path, *args, **kwargs) -> int`

Save the field directly without endianness or precision conversion.

#### `load_direct(path, *args, **kwargs) -> int`

Load the field directly. Geometry and multiplicity are determined during
loading. Returns the number of bytes read (0 on failure).

#### `save_64(path, *args, **kwargs) -> int`

Save with 64-bit endianness conversion:
- Local file → big-endian 64-bit
- Shuffled archive → little-endian 64-bit

#### `load_64(path, *args, **kwargs) -> int`

Load with 64-bit endianness conversion (inverse of `save_64`).

#### `save_double(path, *args, **kwargs) -> int`

Alias for `save_64`.

#### `load_double(path, *args, **kwargs) -> int`

Alias for `load_64`.

#### `save_float_from_double(path, *args, **kwargs) -> int`

Convert a double-precision field to single-precision, apply endianness
conversion, and save.

#### `load_double_from_float(path, *args, **kwargs) -> int`

Load a single-precision field and convert to double-precision.

### Endianness and Precision Conversion

#### `to_from_endianness(tag)`

Convert between native endianness and the specified endianness in-place.

| Tag | Meaning |
|---|---|
| `"big_32"` | Big-endian, 32-bit elements |
| `"big_64"` | Big-endian, 64-bit elements |
| `"little_32"` | Little-endian, 32-bit elements |
| `"little_64"` | Little-endian, 64-bit elements |

#### `float_from_double(f: FieldBase)`

Convert double-precision field `f` to single-precision in `self` (which must
be `FieldRealF`).

#### `double_from_float(ff: FieldRealF)`

Convert single-precision field `ff` to double-precision in `self`.

### Checksums

#### `crc32() -> int`

Compute a CRC-32 checksum of the raw field data.

### Pickle Support

`FieldBase` supports `pickle` via `__getstate__` / `__setstate__`. This only
works reliably on a single MPI node (or when all nodes hold identical data).

---

## `SelectedFieldBase` Class

`SelectedFieldBase` is the abstract base class for all concrete
`SelectedField` types. It mirrors `FieldBase` but operates only on sites
selected by a `FieldSelection`.

Key differences from `FieldBase`:

- Initialized with a `FieldSelection` (`fsel`) instead of a `Geometry`.
- `__len__` returns `n_elems` (number of selected sites), not total sites.
- Arithmetic operators (`+=`, `-=`, `*=`) work only with other
  `SelectedFieldBase` instances.
- Provides `glb_sum_tslice(t_dir=3)` — globally sum the selected field over
  time slices, returning a `SelectedPoints`.

Methods `cast_from`, `get_data_sig`, `__getitem__`/`__setitem__`,
`get_elems`/`set_elems`, `get_elem`/`set_elem`, all `save_*`/`load_*`,
`to_from_endianness`, `float_from_double`, `double_from_float`, and pickle
support all behave analogously to `FieldBase` (see above).

---

## `SelectedPointsBase` Class

`SelectedPointsBase` is the abstract base class for all concrete
`SelectedPoints` types. It holds data at a discrete set of lattice points
selected by a `PointsSelection`.

Key differences from `FieldBase`:

- Initialized with a `PointsSelection` (`psel`) instead of a `Geometry`.
- `__len__` returns `n_points`.
- Provides `save_str()` / `load_str(content)` for text-based serialization.
- Provides `to_numpy()` / `from_numpy(arr)` for direct NumPy interop.
- `get_data_sig` normalizes the result by `geo.num_node` when the points
  distribution type is `"g"` (global).
- No arithmetic operators defined.

Methods `cast_from`, `get_data_sig`, `__getitem__`/`__setitem__`,
`get_elems`/`set_elems`, `get_elem`/`set_elem`, and pickle support all behave
analogously to `FieldBase` (see above).

---

## Field Manipulation Utilities

### `split_fields`

```python
split_fields(fs, f)
```

Split field `f` into a list of fields `fs` by the last axis of the
multiplicity. Each element of `fs` can be a pre-existing `FieldBase` (must
have matching `ctype`) or `None` (in which case a new `Field` is created).

### `merge_fields`

```python
merge_fields(f, fs)
```

Merge a list of fields `fs` into field `f` by concatenating along the
multiplicity axis. Inverse of `split_fields`.

### `merge_fields_ms`

```python
merge_fields_ms(f, fms)
```

Merge selected multiplicity components. `fms` is a list of `(field, m)` pairs.
For each output multiplicity index `m`, the data comes from `fms[m][0]` at
multiplicity `fms[m][1]`:

```python
f.get_elem(x, m) == fms[m][0].get_elem(x, fms[m][1])
```

### `mk_merged_fields_ms`

```python
mk_merged_fields_ms(fms) -> FieldBase
```

Convenience wrapper that creates and returns the merged field. Equivalent to:

```python
f = Field(ctype)
merge_fields_ms(f, fms)
return f
```

---

## Examples

### Creating and Populating a Field

```python
import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
f = q.Field(q.ElemTypeRealD, geo, 4)

# Fill with random data
rng = q.RngState("seed")
f.set_rand(rng, 1.0, -1.0)

# Access as NumPy array
arr = np.asarray(f)
print(arr.shape)  # (n_sites, 4)

q.end_with_mpi()
```

### Arithmetic Operations

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
f1 = q.Field(q.ElemTypeRealD, geo, 1)
f2 = q.Field(q.ElemTypeRealD, geo, 1)

rng = q.RngState("seed")
f1.set_rand(rng, 1.0, -1.0)
f2.set_rand(rng, 1.0, -1.0)

f1 += f2        # element-wise addition
f1 *= 2.0       # scalar multiplication
f1 -= f2

q.end_with_mpi()
```

### Checksumming with `get_data_sig`

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
f = q.Field(q.ElemTypeComplexD, geo, 4)
rng = q.RngState("init")
f.set_rand(rng, 1.0, -1.0)

sig = f.get_data_sig(q.RngState("sig"))
print(f"Signature: {sig}")

q.end_with_mpi()
```

### Type Casting Between Fields

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
f_real = q.Field(q.ElemTypeRealD, geo, 2)
rng = q.RngState("seed")
f_real.set_rand(rng, 1.0, -1.0)

f_complex = q.Field(q.ElemTypeComplexD)
f_complex.cast_from(f_real)  # RealD(mult=2) -> ComplexD(mult=1)

q.end_with_mpi()
```

### Save and Load

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
f = q.Field(q.ElemTypeRealD, geo, 1)
rng = q.RngState("seed")
f.set_rand(rng, 1.0, -1.0)

# Save to file (with endianness conversion)
f.save_double("field.dat")

# Load back
g = q.Field(q.ElemTypeRealD)
g.load_double("field.dat")

q.end_with_mpi()
```

### Splitting and Merging Fields

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
f = q.Field(q.ElemTypeRealD, geo, 6)
rng = q.RngState("seed")
f.set_rand(rng, 1.0, -1.0)

# Split into 2 fields of multiplicity 3
fs = [None, None]
q.split_fields(fs, f)
print(fs[0].multiplicity)  # 3
print(fs[1].multiplicity)  # 3

# Merge back
g = q.Field(q.ElemTypeRealD)
q.merge_fields(g, fs)

q.end_with_mpi()
```

### Merging Selected Components

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
f0 = q.Field(q.ElemTypeRealD, geo, 3)
f1 = q.Field(q.ElemTypeRealD, geo, 3)
rng = q.RngState("seed")
f0.set_rand(rng, 1.0, -1.0)
f1.set_rand(rng, 1.0, -1.0)

# Merge: output[m] comes from f0[0] for m=0, f1[1] for m=1
fms = [(f0, 0), (f1, 1)]
f = q.mk_merged_fields_ms(fms)
print(f.multiplicity)  # 2

q.end_with_mpi()
```
