# `qlat.field_double` — Real-Double Field Operations

Source: `qlat/qlat/field_double.py`

> **Note:** Update this document when updating the source file.

## Outline

- `set_checkers(field)` — set checkerboard pattern on a real-double field.
- `set_double_from_complex(field, cf)` — real part of complex field.
- `set_complex_from_double(field, sf)` — promote real field to complex.
- `set_abs_from_complex(field, cf)` — element-wise absolute value.
- `set_ratio_double(field, sf1, sf2)` — element-wise ratio `sf1 / sf2`.
- `less_than_double(field, sf2, mask)` — element-wise comparison.
- `invert_double(field)` — element-wise inversion `1 / field`.
- `multiply_double(field, factor)` — element-wise field × field product.

## Overview

This module provides element-wise manipulation routines for fields whose
element type is `ElemTypeRealD` (real double precision).  Each function
delegates to the corresponding C++ implementation in `qlat.c` after
verifying element-type compatibility via `assert` statements.

These utilities are typically used for scalar observable fields such as
plaquette densities, topological charge densities, or mask/weight fields.

## API Reference

### `set_checkers(field)`

Fill `field` with a checkerboard pattern (±1 on alternating sites).
Marks as possibly unnecessary in the source.

**Raises:** `AssertionError` if `field.ctype` is not `ElemTypeRealD`.

---

### `set_double_from_complex(field, cf)`

Store the real part of complex field `cf` into real-double field `field`.

| Parameter | Type | Description |
|---|---|---|
| `field` | `FieldBase` | Destination, must be `ElemTypeRealD`. |
| `cf` | `FieldBase` | Source, must be `ElemTypeComplexD`. |

---

### `set_complex_from_double(field, sf)`

Promote real-double field `sf` into complex field `field` (imaginary
part set to zero).

| Parameter | Type | Description |
|---|---|---|
| `field` | `FieldBase` | Destination, must be `ElemTypeComplexD`. |
| `sf` | `FieldBase` | Source, must be `ElemTypeRealD`. |

---

### `set_abs_from_complex(field, cf)`

Store the element-wise absolute value of complex field `cf` into
real-double field `field`.

| Parameter | Type | Description |
|---|---|---|
| `field` | `FieldBase` | Destination, must be `ElemTypeRealD`. |
| `cf` | `FieldBase` | Source, must be `ElemTypeComplexD`. |

---

### `set_ratio_double(field, sf1, sf2)`

Compute the element-wise ratio `sf1 / sf2` and store the result in
`field`.  All three fields must be `ElemTypeRealD`.

| Parameter | Type | Description |
|---|---|---|
| `field` | `FieldBase` | Destination. |
| `sf1` | `FieldBase` | Numerator. |
| `sf2` | `FieldBase` | Denominator. |

---

### `less_than_double(field, sf2, mask)`

For each element, set the corresponding entry in `mask` to 1.0 if
`field[i] < sf2[i]`, otherwise 0.0.  All fields must be
`ElemTypeRealD`.

| Parameter | Type | Description |
|---|---|---|
| `field` | `FieldBase` | Left-hand side of comparison. |
| `sf2` | `FieldBase` | Right-hand side of comparison. |
| `mask` | `FieldBase` | Output mask field. |

---

### `invert_double(field)`

Replace every element of `field` with its reciprocal `1 / element`.

**Raises:** `AssertionError` if `field.ctype` is not `ElemTypeRealD`.

---

### `multiply_double(field, factor)`

Element-wise product `field[i] *= factor[i]`.  Both fields must be
`ElemTypeRealD`.

| Parameter | Type | Description |
|---|---|---|
| `field` | `FieldBase` | Field to be modified in-place. |
| `factor` | `FieldBase` | Multiplicative factor field. |

## Examples

```python
import qlat as q
import qlat_utils as qu

q.begin_with_mpi([1, 1, 1, 4])

geo = q.Geometry([4, 4, 4, 8], 1)
f_real = q.Field(qu.ElemTypeRealD, geo, 1)
f_complex = q.Field(qu.ElemTypeComplexD, geo, 1)
f_real.set_zero()
f_complex.set_zero()

# Convert complex field to real (takes real part)
q.set_double_from_complex(f_real, f_complex)

# Invert all elements in-place
q.invert_double(f_real)

q.end_with_mpi()
```
