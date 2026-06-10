# `qlat_utils.load_prop` — Propagator Loader for Lattice QCD Measurements

Source: `qlat-utils/qlat_utils/load_prop.py`

> **Note:** Update this document when updating the source file.

## Outline

- `load_prop` — convert raw propagator data into `WilsonMatrix` form
- Supports plain `WilsonMatrix` objects, `0`, and AMA-wrapped data
- Handles the g5-hermitian conjugate convention

## Overview

In domain-wall fermion (DWF) lattice QCD calculations, propagators are
often stored in a compressed or intermediate representation and must be
converted into the 12×12 `WilsonMatrix` format before contraction.
Additionally, many measurements require applying the g5-hermitian
conjugate to the propagator.

`load_prop` provides a single entry point for this conversion.  It
transparently handles both plain data and `AmaVal`-wrapped data (via
`ama_apply1`), so the same function works for sloppy, fine, and
multi-accuracy propagators.

## Detailed Sections

### `load_prop(x)`

Convert propagator data `x` into `WilsonMatrix` form.

**Parameters:**

- `x` — raw propagator data.  Can be:
  - A `WilsonMatrix` object, or `0` (accepted by `as_wilson_matrix`)
  - An `AmaVal` wrapping a `WilsonMatrix`
  - A `tuple` `("g5_herm", data)` indicating that the g5-hermitian
    conjugate should be applied; `data` may itself be an `AmaVal`

**Returns:**

- A `WilsonMatrix` (plain value case), or an `AmaVal` containing
  `WilsonMatrix` corrections (AMA case).

**Convention:**

When the input is a 2-tuple `("g5_herm", raw_data)`, the function
applies `as_wilson_matrix_g5_herm` to the raw data, which computes:

```
WM_out = g5 * WM^dagger * g5
```

This is the standard g5-hermiticity relation for quark propagators in
the twisted-mass and domain-wall formulations.

### Dependencies

- `as_wilson_matrix` from `qlat_utils.c` — converts input to a
  `WilsonMatrix` (accepts `WilsonMatrix` or `0`)
- `as_wilson_matrix_g5_herm` from `qlat_utils.c` — converts and applies
  the g5-hermitian conjugate in one step
- `ama_apply1` from `qlat_utils.ama` — applies a unary function through
  all AMA correction terms

## Examples

### Loading a plain propagator

```python
import qlat_utils as q
import numpy as np

# Create a zero WilsonMatrix and load it
raw = q.WilsonMatrix()
raw.set_zero()
wm = q.load_prop(raw)
assert isinstance(wm, q.WilsonMatrix)
# wm is a WilsonMatrix

# Alternatively, pass 0 directly:
wm2 = q.load_prop(0)
assert isinstance(wm2, q.WilsonMatrix)
```

### Loading a g5-hermitian propagator

```python
import qlat_utils as q
import numpy as np

raw = q.WilsonMatrix()
raw.set_zero()
wm = q.load_prop(("g5_herm", raw))
assert isinstance(wm, q.WilsonMatrix)
# wm = g5 * raw^dagger * g5  as a WilsonMatrix
```

### Loading an AMA-wrapped propagator

```python
import qlat_utils as q
import numpy as np

# Create WilsonMatrix objects to wrap in AMA
wm_sloppy = q.WilsonMatrix()
wm_sloppy.set_zero()
wm_fine = q.WilsonMatrix()
wm_fine.set_zero()
np.asarray(wm_fine)[:] = 0.01 * np.ones((12, 12), dtype=np.complex128)

v = q.mk_ama_val(
    wm_sloppy,
    ("point", (0, 0, 0, 0)),
    [wm_sloppy, wm_fine],
    [0, 1],
    [1.0, 0.1],
)
wm_ama = q.load_prop(v)
assert isinstance(wm_ama, q.AmaVal)
# wm_ama is an AmaVal containing WilsonMatrix corrections
result = q.ama_extract(wm_ama)
assert isinstance(result, q.WilsonMatrix)
```

### Loading a g5-hermitian AMA propagator

```python
import qlat_utils as q
import numpy as np

wm_sloppy = q.WilsonMatrix()
wm_sloppy.set_zero()
wm_fine = q.WilsonMatrix()
wm_fine.set_zero()
np.asarray(wm_fine)[:] = 0.01 * np.ones((12, 12), dtype=np.complex128)

v = q.mk_ama_val(
    wm_sloppy,
    ("point", (0, 0, 0, 0)),
    [wm_sloppy, wm_fine],
    [0, 1],
    [1.0, 0.1],
)
wm_ama = q.load_prop(("g5_herm", v))
assert isinstance(wm_ama, q.AmaVal)
# Each AMA correction is independently g5-hermitian conjugated
result = q.ama_extract(wm_ama)
assert isinstance(result, q.WilsonMatrix)
```
