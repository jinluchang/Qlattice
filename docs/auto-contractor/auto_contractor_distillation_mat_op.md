# `auto_contractor.distillation_mat_op` — Distillation Matrix Operations

Source: `qlat/auto_contractor/distillation_mat_op.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Einsum Cache Path](#einsum-cache-path)
3. [Gamma Matrix](#gamma-matrix)
4. [Wilson Matrix Operations](#wilson-matrix-operations)
5. [Spin Matrix Operations](#spin-matrix-operations)
6. [Colour Matrix Operations](#colour-matrix-operations)
7. [Cross-Type Trace Operations](#cross-type-trace-operations)
8. [Cross-Type Multiplication Operations](#cross-type-multiplication-operations)
9. [Examples](#examples)

---

## Overview

`distillation_mat_op` provides NumPy-based matrix operations for the
distillation framework, where propagators are stored as dense NumPy arrays
indexed by `(spin, noise_vector)` rather than as compact `WilsonMatrix` objects.
All operations use `np.einsum` with cached contraction paths for performance.

The module distinguishes three matrix types by their index structure:

| Abbreviation | Type | Index Convention |
|---|---|---|
| `wm` | Wilson matrix | `(spin_snk, noise_snk, spin_src, noise_src)` — 4-index |
| `sm` | Spin matrix | `(spin_snk, spin_src)` — 2-index |
| `cm` | Colour/noise matrix | `(noise_snk, noise_src)` — 2-index |

## Einsum Cache Path

### `einsum_cache_path(subscripts, *operands, optimize=None)`

Wrapper around `np.einsum` that supports cached contraction paths.  When
`optimize` is an empty list `[]`, the optimal path is computed via
`np.einsum_path` with `optimize="optimal"`, stored into the list, and printed.
Subsequent calls reuse the cached path.

## Gamma Matrix

### `get_gamma_matrix(mu)`

Return the gamma matrix for direction `mu` (0–3 for spatial, 5 for γ₅) as a
NumPy array.  Results are cached via `@functools.lru_cache`.

### `wilson_matrix_g5_herm(x)`

Compute `γ₅ · x† · γ₅` for a Wilson matrix `x`.  Used to obtain the
G-parity-conjugate propagator.

### `load_prop(x)`

Load a propagator, applying `g5_herm` conjugation if `x` is a tuple
`("g5_herm", prop_data)`.  Uses `q.ama_apply1` for AMA support.

## Wilson Matrix Operations

All Wilson matrix trace/multiply functions use 4-index Einstein summation.

| Function | Einsum | Description |
|---|---|---|
| `mat_tr_wm(mat)` | `iaia->` | Trace of a Wilson matrix |
| `mat_tr_wm_wm(mat1, mat2)` | `iajb,jbia->` | Trace of product of two Wilson matrices |
| `mat_mul_wm_wm(mat1, mat2)` | `iajb,jbkc->iakc` | Product of two Wilson matrices |

## Spin Matrix Operations

| Function | Einsum | Description |
|---|---|---|
| `mat_tr_sm(mat)` | `ii->` | Trace of a spin matrix |
| `mat_tr_sm_sm(mat1, mat2)` | `ij,ji->` | Trace of product of two spin matrices |
| `mat_mul_sm_sm(mat1, mat2)` | `ij,jk->ik` | Product of two spin matrices |

## Colour Matrix Operations

| Function | Einsum | Description |
|---|---|---|
| `mat_tr_cm(mat)` | `ii->` | Trace of a colour matrix |
| `mat_tr_cm_cm(mat1, mat2)` | `ab,ba->` | Trace of product of two colour matrices |
| `mat_mul_cm_cm(mat1, mat2)` | `ab,bc->ac` | Product of two colour matrices |

## Cross-Type Trace Operations

These trace functions combine operators of different index types.

| Function | Einsum | Description |
|---|---|---|
| `mat_tr_wm_sm(mat1, mat2)` | `iaja,ji->` | tr(WilsonMatrix · SpinMatrix) |
| `mat_tr_sm_wm(mat1, mat2)` | — | tr(SpinMatrix · WilsonMatrix), delegates to `mat_tr_wm_sm` |
| `mat_tr_wm_cm(mat1, mat2)` | `iaib,ba->` | tr(WilsonMatrix · ColorMatrix) |
| `mat_tr_cm_wm(mat1, mat2)` | — | tr(ColorMatrix · WilsonMatrix), delegates to `mat_tr_wm_cm` |

## Cross-Type Multiplication Operations

| Function | Einsum | Description |
|---|---|---|
| `mat_mul_wm_sm(mat1, mat2)` | `iajb,jk->iakb` | WilsonMatrix × SpinMatrix |
| `mat_mul_sm_wm(mat1, mat2)` | `ij,jakb->iakb` | SpinMatrix × WilsonMatrix |
| `mat_mul_wm_cm(mat1, mat2)` | `iajb,bc->iajc` | WilsonMatrix × ColorMatrix |
| `mat_mul_cm_wm(mat1, mat2)` | `ab,ibjc->iajc` | ColorMatrix × WilsonMatrix |

Each cross-type function has a corresponding `einsum_optimize_*` list that
caches the contraction path on first use.

## Examples

```python
import qlat as q
import numpy as np

q.begin_with_mpi([[1, 1, 1, 4]])

from auto_contractor.distillation_mat_op import (
    get_gamma_matrix,
    mat_mul_wm_wm,
    mat_tr_wm,
    einsum_cache_path,
)

# Create random 4-index Wilson matrices (spin=4, noise=10)
rng = np.random.default_rng(42)
ns, nc = 4, 10
wm1 = rng.standard_normal((ns, nc, ns, nc)) + 1j * rng.standard_normal((ns, nc, ns, nc))
wm2 = rng.standard_normal((ns, nc, ns, nc)) + 1j * rng.standard_normal((ns, nc, ns, nc))

# Multiply and trace
wm3 = mat_mul_wm_wm(wm1, wm2)
tr = mat_tr_wm(wm3)
print(f"tr(S1 * S2) = {tr}")

# Gamma matrix
g5 = get_gamma_matrix(5)
print(f"gamma_5 shape: {g5.shape}")

q.end_with_mpi()
```
