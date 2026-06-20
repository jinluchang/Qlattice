# `auto_contractor.runtime` — Runtime for Compiled Expression Evaluation

Source: `qlat/auto_contractor/runtime.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Exported Names](#exported-names)
   - [NumPy and Timers](#numpy-and-timers)
   - [AMA Utilities](#ama-utilities)
   - [Matrix Types](#matrix-types)
   - [Matrix Operations](#matrix-operations)
   - [Gamma and Propagator Utilities](#gamma-and-propagator-utilities)
3. [Usage](#usage)
4. [Examples](#examples)

---

## Overview

`runtime` re-exports the symbols needed at runtime by compiled auto-contractor
expressions (`CExpr`). When a `CExpr` is evaluated, it uses functions and types
from this module to load propagators, perform matrix arithmetic, and compute
traces. Importing this module ensures all required runtime dependencies are
available.

## Exported Names

### NumPy and Timers

| Name | Source | Description |
|---|---|---|
| `np` | `numpy` | NumPy array library |
| `timer` | `qlat.timer` | Scope-based timer decorator |
| `timer_flops` | `qlat.timer_flops` | Timer with FLOPS counting |

### AMA Utilities

| Name | Source | Description |
|---|---|---|
| `ama_list` | `qlat_utils.ama` | Accumulate multi-shift solver results |
| `ama_apply1` | `qlat_utils.ama` | Apply AMA correction to a single result |
| `ama_counts` | `qlat_utils.ama` | Count AMA samples |
| `ama_extract` | `qlat_utils.ama` | Extract corrected value from AMA list |

### Matrix Types

| Name | Source | Description |
|---|---|---|
| `WilsonMatrix` | `qlat_utils` | $4\times3$ spin-color matrix (propagator element) |
| `SpinMatrix` | `qlat_utils` | $4\times4$ spin matrix (gamma matrices) |
| `ColorMatrix` | `qlat_utils` | $3\times3$ color matrix (gauge links) |

### Matrix Operations

All `mat_*` functions are low-level C++ routines exposed via `qlat_utils`:

**Trace operations:**

| Function | Description |
|---|---|
| `mat_tr_sm(m)` | Trace of SpinMatrix |
| `mat_tr_cm(m)` | Trace of ColorMatrix |
| `mat_tr_wm(m)` | Trace of WilsonMatrix |
| `mat_tr_wm_wm(a, b)` | Trace of product $\mathrm{tr}(A \cdot B)$, both WilsonMatrix |
| `mat_tr_wm_sm(a, b)` | Trace of WilsonMatrix times SpinMatrix |
| `mat_tr_sm_wm(a, b)` | Trace of SpinMatrix times WilsonMatrix |
| `mat_tr_sm_sm(a, b)` | Trace of SpinMatrix times SpinMatrix |
| `mat_tr_wm_cm(a, b)` | Trace of WilsonMatrix times ColorMatrix |
| `mat_tr_cm_wm(a, b)` | Trace of ColorMatrix times WilsonMatrix |
| `mat_tr_cm_cm(a, b)` | Trace of ColorMatrix times ColorMatrix |

**Multiplication operations:**

| Function | Description |
|---|---|
| `mat_mul_a_wm(a, m)` | Scalar times WilsonMatrix |
| `mat_mul_a_sm(a, m)` | Scalar times SpinMatrix |
| `mat_mul_a_cm(a, m)` | Scalar times ColorMatrix |
| `mat_mul_wm_wm(a, b)` | WilsonMatrix times WilsonMatrix |
| `mat_mul_wm_sm(a, b)` | WilsonMatrix times SpinMatrix |
| `mat_mul_sm_wm(a, b)` | SpinMatrix times WilsonMatrix |
| `mat_mul_sm_sm(a, b)` | SpinMatrix times SpinMatrix |
| `mat_mul_wm_cm(a, b)` | WilsonMatrix times ColorMatrix |
| `mat_mul_cm_wm(a, b)` | ColorMatrix times WilsonMatrix |
| `mat_mul_cm_cm(a, b)` | ColorMatrix times ColorMatrix |

**Addition operations:**

| Function | Description |
|---|---|
| `mat_add_wm_wm(a, b)` | WilsonMatrix + WilsonMatrix |
| `mat_add_sm_sm(a, b)` | SpinMatrix + SpinMatrix |
| `mat_add_cm_cm(a, b)` | ColorMatrix + ColorMatrix |

**Special operations:**

| Function | Description |
|---|---|
| `mat_epsilon_contraction_wm_wm_wm(a, b, c)` | $\epsilon_{ijk}$ contraction of three WilsonMatrix objects |

### Gamma and Propagator Utilities

| Function | Description |
|---|---|
| `load_prop(...)` | Load a propagator from disk |
| `get_gamma_matrix(tag)` | Get gamma matrix by tag (0–3, 5) |
| `wilson_matrix_g5_herm(m)` | Compute $\gamma_5 m^\dagger \gamma_5$ |

### Auto-contractor Functions

| Name | Description |
|---|---|
| `aff` | `auto_contractor.auto_fac_funcs` module — callable functions for `CExpr` evaluation |

## Usage

The runtime module is typically imported by compiled expression code, not
directly by users:

```python
from qlat.auto_contractor.runtime import *
```

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])

from qlat.auto_contractor.runtime import (
    WilsonMatrix, SpinMatrix, get_gamma_matrix,
    mat_tr_wm_wm, mat_mul_wm_sm,
)

# Create a WilsonMatrix and multiply by gamma_5
wm = WilsonMatrix()
g5 = get_gamma_matrix(5)
result = mat_mul_wm_sm(wm, g5)
trace = mat_tr_wm_wm(result, wm)

q.end_with_mpi()
```
