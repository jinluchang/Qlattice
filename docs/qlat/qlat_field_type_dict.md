# `qlat.field_type_dict` — Field Element-Type Registry

Source: `qlat/qlat/field_type_dict.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Type Dictionaries](#type-dictionaries)
   - [`field_type_dict`](#field_type_dict)
   - [`selected_field_type_dict`](#selected_field_type_dict)
   - [`selected_points_type_dict`](#selected_points_type_dict)
3. [Element-Type Category Lists](#element-type-category-lists)
   - [`field_ctypes_complex`](#field_ctypes_complex)
   - [`field_ctypes_double`](#field_ctypes_double)
   - [`field_ctypes_float`](#field_ctypes_float)
   - [`field_ctypes_long`](#field_ctypes_long)
   - [`field_ctypes_char`](#field_ctypes_char)
4. [Internal Lists](#internal-lists)
5. [Examples](#examples)

---

## Overview

`field_type_dict` provides a centralised registry of field element types used
by Qlattice's field I/O, serialisation, and dispatch machinery. Other modules
consult these dictionaries and lists to determine which C++ code path to
invoke for a given element type (e.g. `ColorMatrix`, `WilsonVector`,
`ComplexD`, `RealD`, `Long`, `Char`).

The three dictionaries are initially empty and are populated at import time
by the C/C++ extension modules. The category lists group the `ElemType*`
classes by their numeric representation (complex, double, float, long, char)
for convenient lookup.

---

## Type Dictionaries

### `field_type_dict`

```python
field_type_dict: dict
```

Maps `Field` subclass names (or element-type identifiers) to their
corresponding C++ field types. Used by `Field` constructors and I/O routines
to select the correct template instantiation.

### `selected_field_type_dict`

```python
selected_field_type_dict: dict
```

Maps `SelectedField` subclass names to their C++ types. Analogous to
`field_type_dict` but for fields that operate on a subset of lattice sites.

### `selected_points_type_dict`

```python
selected_points_type_dict: dict
```

Maps `SelectedPoints` subclass names to their C++ types. Used for
point-based field data that has been extracted at specific lattice
coordinates.

---

## Element-Type Category Lists

Each list groups `ElemType*` classes from `qlat_utils` by their underlying
numeric representation. These are used internally to dispatch to the
appropriate C++ template specialisation.

### `field_ctypes_complex`

```python
field_ctypes_complex = [
    ElemTypeColorMatrix,
    ElemTypeWilsonMatrix,
    ElemTypeNonRelWilsonMatrix,
    ElemTypeIsospinMatrix,
    ElemTypeSpinMatrix,
    ElemTypeWilsonVector,
    ElemTypeComplexD,
]
```

Element types backed by complex-valued (double-precision) data. Includes all
spin-color matrix types and scalar complex.

### `field_ctypes_double`

```python
field_ctypes_double = [
    ElemTypeColorMatrix,
    ElemTypeWilsonMatrix,
    ElemTypeNonRelWilsonMatrix,
    ElemTypeIsospinMatrix,
    ElemTypeSpinMatrix,
    ElemTypeWilsonVector,
    ElemTypeComplexD,
    ElemTypeRealD,
]
```

Element types stored as double-precision real numbers. Superset of
`field_ctypes_complex` plus `ElemTypeRealD`.

### `field_ctypes_float`

```python
field_ctypes_float = [
    ElemTypeComplexF,
    ElemTypeRealF,
]
```

Single-precision element types.

### `field_ctypes_long`

```python
field_ctypes_long = [
    ElemTypeLong,
    ElemTypeInt64t,
]
```

Integer (64-bit) element types.

### `field_ctypes_char`

```python
field_ctypes_char = [
    ElemTypeChar,
    ElemTypeInt8t,
]
```

Byte-level element types.

---

## Internal Lists

These are not exported in `__all__` but are available within the module:

| Name | Contents |
|---|---|
| `field_ctypes_complex_f` | `[ElemTypeComplexF]` |
| `field_ctypes_complex_or_complex_f` | Union of complex and single-precision complex |
| `field_ctypes_double_or_float` | Union of double and float categories |

---

## Examples

### Inspecting the Type Dictionaries

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

print("field_type_dict keys:")
for key in sorted(q.field_type_dict.keys()):
    print(f"  {key}")

print("\nselected_field_type_dict keys:")
for key in sorted(q.selected_field_type_dict.keys()):
    print(f"  {key}")

print("\nselected_points_type_dict keys:")
for key in sorted(q.selected_points_type_dict.keys()):
    print(f"  {key}")

q.end_with_mpi()
```

### Checking Element-Type Categories

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

print("Complex element types:")
for et in q.field_ctypes_complex:
    print(f"  {et}")

print("\nDouble element types:")
for et in q.field_ctypes_double:
    print(f"  {et}")

print("\nFloat element types:")
for et in q.field_ctypes_float:
    print(f"  {et}")

print("\nLong element types:")
for et in q.field_ctypes_long:
    print(f"  {et}")

print("\nChar element types:")
for et in q.field_ctypes_char:
    print(f"  {et}")

q.end_with_mpi()
```

### Using Type Lists for Dispatch

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

def is_complex_type(elem_type):
    return elem_type in q.field_ctypes_complex

def is_double_type(elem_type):
    return elem_type in q.field_ctypes_double

for et in q.field_ctypes_double:
    print(f"{et}: complex={is_complex_type(et)}, double={is_double_type(et)}")

q.end_with_mpi()
```
