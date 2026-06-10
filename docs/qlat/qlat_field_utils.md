# `qlat.field_utils` — Field Utility Functions

Source: `qlat/qlat/field_utils.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Field Expansion](#field-expansion)
3. [Communication Plans](#communication-plans)
4. [Fast Fourier Transform](#fast-fourier-transform)
5. [Norm and Square Root](#norm-and-square-root)
6. [Coordinate Shifting](#coordinate-shifting)
7. [Field Shuffling](#field-shuffling)
8. [Examples](#examples)

---

## Overview

`field_utils` provides standalone utility functions that operate on lattice
fields. These functions handle halo expansion, communication planning, FFT,
norms, element-wise operations, coordinate shifts, and data redistribution
(shuffling) across MPI nodes.

## Field Expansion

```python
f_e = field_expanded(f, expansion_left, expansion_right)
```

Create a copy of field `f` on an expanded geometry with the given halo
widths on each side. The returned field has the same data on the original
sites and undefined values in the halo.

| Function | Description |
|---|---|
| `field_expanded(f, expansion_left, expansion_right)` | Return a new field on the expanded geometry |
| `refresh_expanded(field, comm_plan=None)` | Refresh the halo data via MPI communication |
| `refresh_expanded_1(field)` | Refresh halo using a single-direction communication pattern |

## Communication Plans

The `CommMarks` and `CommPlan` classes manage field halo communication.

```python
cp = make_field_expand_comm_plan(comm_marks)
```

| Class / Function | Description |
|---|---|
| `CommMarks` | A `Field(ElemTypeInt8t)` marking which sites need communication |
| `CommPlan` | A precomputed communication plan derived from `CommMarks` |
| `make_field_expand_comm_plan(comm_marks)` | Build a `CommPlan` from a `CommMarks` field |

`CommPlan` supports copy and assignment (`@=`).

## Fast Fourier Transform

The `FastFourierTransform` class applies FFT to one or more fields.

```python
fft = mk_fft(is_forward, is_only_spatial=False, is_normalizing=False, mode_fft=1)
fields_out = fft * fields_in   # apply FFT (returns copies)
```

| Class / Function | Description |
|---|---|
| `FastFourierTransform(fft_infos, ...)` | Construct an FFT operator for specified directions |
| `mk_fft(is_forward, ...)` | Factory: create a standard 3D or 4D FFT operator |

Parameters for `mk_fft`:

- `is_forward` — `True` for forward FFT, `False` for inverse.
- `is_only_spatial` — If `True`, transform only directions 0, 1, 2.
- `is_normalizing` — If `True`, apply `1/sqrt(N)` normalization.
- `mode_fft` — FFT implementation mode (0 or 1).

When applied via `*`, `FastFourierTransform` accepts a single `FieldBase` or
a list of `FieldBase` objects. Returns a list (or single field if one was given).

## Norm and Square Root

| Function | Description |
|---|---|
| `qnorm_field(f)` | Per-site squared norm; returns `SelectedFieldRealD` or `SelectedPointsRealD` |
| `sqrt_field(f)` | Element-wise square root; dispatches by field type |

`sqrt_field` supports `FieldRealD`, `SelectedFieldRealD`, and
`SelectedPointsRealD`. The underlying type-specific helpers are
`sqrt_field_real_d`, `sqrt_selected_field_real_d`, and
`sqrt_selected_points_real_d`.

## Coordinate Shifting

```python
sf = field_shift(f, shift)        # shift any FieldBase
sf = field_char_shift(f, shift)   # shift a FieldChar directly
```

| Function | Description |
|---|---|
| `field_shift(f, shift)` | Return a new shifted copy of `f`; `self` is unchanged |
| `field_char_shift(f, shift)` | Low-level shift for `FieldChar` |

`shift` is a `Coordinate` specifying the cyclic displacement. Roughly:
`shifted_field[(xg + shift) % total_site] == original_field[xg]`.

## Field Shuffling

Shuffling redistributes field data across MPI nodes according to a new node
geometry.

```python
f_list = shuffle_field(f, new_size_node)
shuffle_field_back(f, f_list, new_size_node)   # reassemble in-place
```

| Function | Description |
|---|---|
| `shuffle_field(f, new_size_node)` | Split `f` into sub-fields for the new node layout |
| `shuffle_field_back(f, f_list, new_size_node)` | Reassemble `f` from shuffled sub-fields (modifies `f` in-place) |
| `shuffle_field_char(f, new_size_node)` | Low-level shuffle for `FieldChar` |
| `shuffle_field_char_back(fc, f_list, new_size_node)` | Low-level reassemble for `FieldChar` |

## Examples

```python
import qlat as q

q.begin_with_mpi([1, 1, 1, 4])

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))
f = q.FieldRealD(geo, 1)
f.set_zero()

# Expand the field halo by 1 site on each side
f_e = q.field_expanded(f, 1, 1)
print(f"Expanded n_sites: {f_e.n_sites}")

# FFT (forward, spatial only)
fft = q.mk_fft(True, is_only_spatial=True)
f_fft = fft * f

# Norm
norm = f.qnorm()
print(f"qnorm: {norm}")

# Square root
f_sqrt = q.sqrt_field(f)

# Shift by (1,0,0,0)
f_shifted = q.field_shift(f, q.Coordinate([1, 0, 0, 0]))

q.end_with_mpi()
```
