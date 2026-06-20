# `auto_contractor.runtime_distillation` — Runtime for Distillation Evaluation

Source: `qlat/auto_contractor/runtime_distillation.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Differences from `runtime`](#differences-from-runtime)
3. [Exported Names](#exported-names)
4. [Usage](#usage)
5. [Examples](#examples)

---

## Overview

`runtime_distillation` provides the runtime environment for evaluating
compiled auto-contractor expressions using the distillation method. It
exports the same interface as `runtime` but replaces the propagator and
matrix operation backends with distillation-aware implementations from
`auto_contractor.distillation_mat_op`.

Distillation is a technique where quark propagators are projected onto a
subspace of low-mode eigenvectors of the Laplacian, reducing noise in
correlation functions.

## Differences from `runtime`

The key difference is the source of matrix operations:

| Aspect | `runtime` | `runtime_distillation` |
|---|---|---|
| Matrix ops source | `qlat_utils` (C++) | `auto_contractor.distillation_mat_op` |
| `WilsonMatrix` | Included | **Not** included |
| `SpinMatrix` | Included | **Not** included |
| `ColorMatrix` | Included | **Not** included |
| `mat_add_*` | Included | **Not** included |
| `mat_epsilon_contraction_*` | Included | **Not** included |

The distillation backend operates on different data structures (per-operator
eigenvector projections rather than point-to-all propagators), so the
WilsonMatrix/SpinMatrix/ColorMatrix types and certain operations are not
applicable.

## Exported Names

### Common with `runtime`

| Name | Source | Description |
|---|---|---|
| `np` | `numpy` | NumPy array library |
| `timer` | `qlat.timer` | Scope-based timer |
| `timer_flops` | `qlat.timer_flops` | Timer with FLOPS counting |
| `ama_list`, `ama_apply1`, `ama_counts`, `ama_extract` | `qlat_utils.ama` | AMA utilities |
| `aff` | `auto_contractor.auto_fac_funcs` | Callable functions for evaluation |

### Distillation Matrix Operations

All from `auto_contractor.distillation_mat_op`:

| Category | Functions |
|---|---|
| Propagator | `load_prop` |
| Gamma | `get_gamma_matrix`, `wilson_matrix_g5_herm` |
| Trace | `mat_tr_sm`, `mat_tr_cm`, `mat_tr_wm`, `mat_tr_wm_wm`, `mat_tr_wm_sm`, `mat_tr_sm_wm`, `mat_tr_sm_sm`, `mat_tr_wm_cm`, `mat_tr_cm_wm`, `mat_tr_cm_cm` |
| Multiply | `mat_mul_wm_wm`, `mat_mul_wm_sm`, `mat_mul_sm_wm`, `mat_mul_sm_sm`, `mat_mul_wm_cm`, `mat_mul_cm_wm`, `mat_mul_cm_cm` |

## Usage

Import this module instead of `runtime` when evaluating compiled expressions
with distillation operators:

```python
from qlat.auto_contractor.runtime_distillation import *
```

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])

from qlat.auto_contractor.runtime_distillation import (
    get_gamma_matrix, mat_tr_wm_wm, load_prop, aff,
)

# Use distillation-aware runtime for expression evaluation
g5 = get_gamma_matrix(5)

q.end_with_mpi()
```
