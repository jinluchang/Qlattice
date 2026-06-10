# `qlat.selected_field_types` — Typed Selected-Field Classes

Source: `qlat/qlat/selected_field_types.pyx.in`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Generated Types](#generated-types)
3. [Construction](#construction)
4. [Properties](#properties)
5. [Assignment and Copying](#assignment-and-copying)
6. [Mutation Operations](#mutation-operations)
7. [Shifting and Reflection](#shifting-and-reflection)
8. [I/O Operations](#io-operations)
9. [NumPy Buffer Protocol](#numpy-buffer-protocol)
10. [Serialization](#serialization)
11. [Selected-Field Type Dictionary](#selected-field-type-dictionary)
12. [Examples](#examples)

---

## Overview

`selected_field_types` is a template module (`.pyx.in`) that generates concrete
`SelectedField{Type}` subclasses of `SelectedFieldBase` for every lattice data
type. Each generated class stores field data only on a subset of lattice sites
identified by a `FieldSelection`, reducing memory and communication costs when
only part of the lattice is needed.

All selected-field classes are registered in `selected_field_type_dict`:
`selected_field_type_dict[ElemTypeRealD]` returns `SelectedFieldRealD`.

## Generated Types

The template generates `SelectedField` variants for all 16 element types:

| Category | Classes |
|---|---|
| Complex matrices | `SelectedFieldColorMatrix`, `SelectedFieldWilsonMatrix`, `SelectedFieldNonRelWilsonMatrix`, `SelectedFieldIsospinMatrix`, `SelectedFieldSpinMatrix` |
| Complex scalars | `SelectedFieldWilsonVector`, `SelectedFieldComplexD`, `SelectedFieldComplexF` |
| Real scalars | `SelectedFieldRealD`, `SelectedFieldRealF` |
| Integer types | `SelectedFieldLong`, `SelectedFieldInt`, `SelectedFieldInt64t`, `SelectedFieldInt32t`, `SelectedFieldInt8t`, `SelectedFieldChar` |

## Construction

```python
sf = SelectedFieldRealD()                   # empty (uninitialized)
sf = SelectedFieldRealD(fsel)               # allocate with field selection
sf = SelectedFieldRealD(fsel, multiplicity) # with given multiplicity
```

`fsel` is a `FieldSelection` object specifying which lattice sites are selected.

## Properties

| Property | Type | Description |
|---|---|---|
| `n_elems` | `int` | Number of selected sites |
| `total_site` | `Coordinate` | Total lattice dimensions |
| `multiplicity` | `int` | Number of elements per site |
| `sizeof_m` | `int` | Size of one element in bytes |
| `geo` | `Geometry` | The lattice geometry |

## Assignment and Copying

The `@=` operator assigns data. It does **not** change `self.fsel`.

```python
sf1 @= sf2                  # SelectedField ← SelectedField (fsel may differ)
sf  @= f                    # SelectedField ← Field (scatter into selected sites)
sf  @= sp                   # SelectedField ← SelectedPoints
```

| Method | Description |
|---|---|
| `copy(is_copying_data=True)` | Deep copy; returns a new `SelectedField` with same `fsel` |
| `set_zero()` | Set all elements to zero |
| `swap(f1)` | Swap contents with another `SelectedField` of the same type <!-- TODO: swap() has a source bug — uses self.psel instead of self.fsel in selected_field_types.pyx.in:142-143 --> |
| `swap_cast(f1)` | Swap raw bytes with a `SelectedFieldChar` |
| `swap_sp_cast(sp, geo)` | Swap raw bytes with a `SelectedPointsChar` (does not affect `psel`/`fsel`) |

## Mutation Operations

```python
sf.set_zero()                                  # zero all elements
sf.set_rand(rng, upper=1.0, lower=0.0)        # uniform random
sf.set_rand_g(rng, center=0.0, sigma=1.0)     # Gaussian random
```

| Method | Description |
|---|---|
| `qnorm()` | Global squared norm (scalar float) |
| `qnorm_field()` | Per-site squared norm → `SelectedFieldRealD` |

## Shifting and Reflection

```python
sf_shifted = sf.shift(shift, is_reflect=False)
```

Return a new shifted `SelectedField`. `self` is not modified. If `shift` is
`None` and `is_reflect` is `False`, returns a plain copy.

## I/O Operations

```python
total_bytes = sf.read_direct(filename)
total_bytes = sf.write_direct(filename)

total_bytes = sf.read_sfr_direct(sfr, filename)          # from ShuffledFieldsReader
total_bytes = sf.write_sfw_direct(sfw, filename)          # via ShuffledFieldsWriter
sf.write_sfw_direct(sfw, filename, skip_if_exist=True)   # skip if already written
```

The shuffled I/O methods (`read_sfr_direct`, `write_sfw_direct`) use a
`ShuffledBitSet` cache from the reader/writer for efficient data transfer.

## NumPy Buffer Protocol

All selected-field types support zero-copy NumPy access:

```python
arr = np.asarray(sf)    # shape: (n_elems, multiplicity, *element_shape)
```

The buffer protocol follows the same conventions as `Field{Type}` but with
`n_elems` (selected sites) instead of `n_sites` (local sites).

## Serialization

```python
state = sf.__getstate__()     # serialize (single-node only)
sf.__setstate__(state)        # deserialize (single-node only)
```

Serialization only works when running on a single node or when all nodes hold
identical data.

## Selected-Field Type Dictionary

```python
from qlat.field_type_dict import selected_field_type_dict

SelectedFieldRealD = selected_field_type_dict[ElemTypeRealD]
```

## Examples

```python
import qlat as q
import numpy as np

q.begin_with_mpi([1, 1, 1, 4])

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))

# Create a field and populate it
f = q.FieldRealD(geo, 1)
rng = q.RngState("seed-1")
f.set_rand_g(rng)

# Build a field selection (select all sites)
fsel = q.FieldSelection(geo, 0)

# Create a selected field
sf = q.SelectedFieldRealD(fsel, 1)
sf @= f                    # scatter from full field
print(f"n_elems: {sf.n_elems}")

# Access as NumPy array
arr = np.asarray(sf)
print(f"Shape: {arr.shape}")

# Norm
norm = sf.qnorm()
print(f"qnorm: {norm}")

# Copy and shift
sf_copy = sf.copy()
sf_shifted = sf.shift(q.Coordinate([1, 0, 0, 0]))

q.end_with_mpi()
```
