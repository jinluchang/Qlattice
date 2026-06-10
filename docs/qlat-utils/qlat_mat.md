# `qlat_utils.mat` — Matrix Types and Operations for Lattice QCD

Source: `qlat-utils/qlat_utils/mat.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Matrix Classes](#matrix-classes)
   - [WilsonMatrix](#wilsonmatrix)
   - [SpinMatrix](#spinmatrix)
   - [ColorMatrix](#colormatrix)
   - [Common Interface](#common-interface)
3. [Trace Functions](#trace-functions)
4. [Multiplication Functions](#multiplication-functions)
5. [Addition Functions](#addition-functions)
6. [Epsilon Contraction](#epsilon-contraction)
7. [Gamma Matrices](#gamma-matrices)
8. [Helper Functions](#helper-functions)
9. [Examples](#examples)

---

## Overview

`mat` provides the three core matrix types used throughout Qlattice for lattice QCD calculations:

- **`WilsonMatrix`** — a 12x12 complex matrix (4 spin x 3 color degrees of freedom) used for quark propagators.
- **`SpinMatrix`** — a 4x4 complex matrix acting on spin indices only.
- **`ColorMatrix`** — a 3x3 complex matrix acting on color indices only (e.g., gauge links).

Each class supports arithmetic operators, NumPy buffer protocol, copying, and conjugation/transpose/adjoint operations. The module also provides free functions for trace, multiplication, addition, gamma matrices, and epsilon contractions.

Python access:

```python
import qlat_utils as q
```

---

## Matrix Classes

### WilsonMatrix

A 12x12 complex double matrix representing a full spin-color propagator.

```python
wm = q.WilsonMatrix()      # zero-initialized
wm = q.as_wilson_matrix(0)  # also zero-initialized
```

### SpinMatrix

A 4x4 complex double matrix acting on spin indices.

```python
sm = q.SpinMatrix()  # zero-initialized
```

### ColorMatrix

A 3x3 complex double matrix acting on color indices (e.g., SU(3) gauge links).

```python
cm = q.ColorMatrix()  # zero-initialized
```

### Common Interface

All three matrix classes share the same interface:

### Copying

| Method | Description |
|---|---|
| `copy(is_copying_data=True)` | Return a copy. If `False`, return a zero matrix. |
| `__copy__()` | Shallow copy (delegates to `copy()`). |
| `__deepcopy__(memo)` | Deep copy (delegates to `copy()`). |
| `__imatmul__(v)` | In-place assignment: `m1 @= m2`. |

### Element Access

| Method | Description |
|---|---|
| `__getitem__(idx)` | Index via `np.asarray(self)[idx]`. |
| `__setitem__(idx, val)` | Assign via `np.asarray(self)[idx] = val`. |

The buffer protocol is implemented, so `np.asarray(wm)` returns a view of shape `(12, 12)` for `WilsonMatrix`, `(4, 4)` for `SpinMatrix`, or `(3, 3)` for `ColorMatrix`.

### Arithmetic

| Operator | Description |
|---|---|
| `+` | Matrix addition |
| `-` | Matrix subtraction |
| `*` (scalar) | Scalar multiplication by `complex` |
| `+=` | In-place addition |
| `-=` | In-place subtraction |
| `*=` | In-place scalar multiplication |
| `==` | Equality check (via `qnorm` of difference) |

### Linear Algebra

| Method | Description |
|---|---|
| `conjugate()` | Element-wise complex conjugate. |
| `transpose()` | Matrix transpose. |
| `adjoint()` | Conjugate transpose (Hermitian adjoint). |
| `g5_herm()` | In-place g5-Hermitization: `M -> g5 * M^adjoint * g5`. (WilsonMatrix only) |
| `qnorm()` | Squared Frobenius norm: `sum |M_ij|^2`. |
| `set_zero()` | Set all elements to zero. |

### Representation

| Method | Description |
|---|---|
| `__repr__()` | Returns a string like `WilsonMatrix([[...]])`. |

---

## Trace Functions

All trace functions return a `complex` value.

| Function | Description |
|---|---|
| `mat_tr_sm(sm)` | `Tr(sm)` |
| `mat_tr_cm(cm)` | `Tr(cm)` |
| `mat_tr_wm(wm)` | `Tr(wm)` |
| `mat_tr_wm_wm(wm1, wm2)` | `Tr(wm1 * wm2)` |
| `mat_tr_wm_sm(wm, sm)` | `Tr(wm * sm)` |
| `mat_tr_sm_wm(sm, wm)` | `Tr(sm * wm)` |
| `mat_tr_sm_sm(sm1, sm2)` | `Tr(sm1 * sm2)` |
| `mat_tr_wm_cm(wm, cm)` | `Tr(wm * cm)` |
| `mat_tr_cm_wm(cm, wm)` | `Tr(cm * wm)` |
| `mat_tr_cm_cm(cm1, cm2)` | `Tr(cm1 * cm2)` |

---

## Multiplication Functions

All multiplication functions return a new matrix.

| Function | Return Type | Description |
|---|---|---|
| `mat_mul_a_wm(a, wm)` | WilsonMatrix | Scalar `a` times `wm` |
| `mat_mul_a_sm(a, sm)` | SpinMatrix | Scalar `a` times `sm` |
| `mat_mul_a_cm(a, cm)` | ColorMatrix | Scalar `a` times `cm` |
| `mat_mul_wm_wm(wm1, wm2)` | WilsonMatrix | `wm1 * wm2` |
| `mat_mul_sm_wm(sm, wm)` | WilsonMatrix | `sm * wm` |
| `mat_mul_wm_sm(wm, sm)` | WilsonMatrix | `wm * sm` |
| `mat_mul_sm_sm(sm1, sm2)` | SpinMatrix | `sm1 * sm2` |
| `mat_mul_cm_wm(cm, wm)` | WilsonMatrix | `cm * wm` |
| `mat_mul_wm_cm(wm, cm)` | WilsonMatrix | `wm * cm` |
| `mat_mul_cm_cm(cm1, cm2)` | ColorMatrix | `cm1 * cm2` |

---

## Addition Functions

| Function | Return Type | Description |
|---|---|---|
| `mat_add_wm_wm(wm1, wm2)` | WilsonMatrix | `wm1 + wm2` |
| `mat_add_sm_sm(sm1, sm2)` | SpinMatrix | `sm1 + sm2` |
| `mat_add_cm_cm(cm1, cm2)` | ColorMatrix | `cm1 + cm2` |

---

## Epsilon Contraction

```python
mat_epsilon_contraction_wm_wm_wm(v_s1, b_s1, v_s2, b_s2, v_s3, b_s3, wm1, wm2, wm3) -> complex
```

Compute the epsilon (Levi-Civita) contraction of three `WilsonMatrix` objects with specified spin index pairs. Parameters `v_s*` are the "vertex" spin indices and `b_s*` are the "baryon" spin indices (each in 0..3). Returns a `complex` scalar.

---

## Gamma Matrices

### `get_gamma_matrix(mu) -> SpinMatrix`

Return the gamma matrix for direction `mu` (0, 1, 2, 3, or 5).

```python
g5 = q.get_gamma_matrix(5)
```

<!-- TODO: gamma_matrix_0..3, gamma_matrix_5 are defined in mat.pyx but declared cdef in mat.pxd (C-level only, not accessible from Python). Consider exposing them as Python-level module constants if user access is desired. -->

---

## Helper Functions

### `as_wilson_matrix(x) -> WilsonMatrix`

Convert `x` to a `WilsonMatrix`. If `x` is already a `WilsonMatrix`, return it directly. If `x == 0`, return a zero `WilsonMatrix`.

### `wilson_matrix_g5_herm(wm) -> WilsonMatrix`

Return `g5 * wm^adjoint * g5` (non-mutating version of `wm.g5_herm()`).

### `as_wilson_matrix_g5_herm(x) -> WilsonMatrix`

Convert `x` to a `WilsonMatrix` and apply `g5_herm`. If `x == 0`, return a zero matrix.

### `benchmark_matrix_functions(count)`

Run a matrix-operation benchmark `count` times. Useful for performance testing.

---

## Examples

### Basic Matrix Operations

```python
import qlat_utils as q

wm1 = q.WilsonMatrix()
wm2 = q.WilsonMatrix()

# Arithmetic
wm3 = wm1 + wm2
wm1 += wm2
wm1 -= wm2
wm3 = q.mat_mul_wm_wm(wm1, wm2)
```

### Working with NumPy

```python
import qlat_utils as q
import numpy as np

sm = q.SpinMatrix()
arr = np.asarray(sm)        # view as (4,4) complex128 array
arr[0, 0] = 1.0 + 0j       # modifies sm in-place
sm[1, 2] = 2.0 + 1j        # also modifies via __setitem__
```

### Trace and Adjoint

```python
import qlat_utils as q

wm = q.WilsonMatrix()
tr = q.mat_tr_wm(wm)
wm_adj = wm.adjoint()
wm.g5_herm()
```

### Gamma Matrices

```python
import qlat_utils as q

g5 = q.get_gamma_matrix(5)
sm = q.SpinMatrix()
g5_sm = q.mat_mul_sm_sm(g5, sm)
```
