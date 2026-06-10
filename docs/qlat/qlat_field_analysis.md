# `qlat.field_analysis` — Field Smearing, Sphere Summation, and Convolution

Source: `qlat/qlat/field_analysis.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Shift Index Utilities](#shift-index-utilities)
   - [`mk_shift_xg_idx_arr`](#mk_shift_xg_idx_arr)
3. [Step-Based Smearing](#step-based-smearing)
   - [`smear_field_step_local`](#smear_field_step_local)
   - [`smear_field_step`](#smear_field_step)
4. [Momentum-Space Kernels](#momentum-space-kernels)
   - [`mk_smear_mom_kernel`](#mk_smear_mom_kernel)
   - [`mk_spatial_smear_mom_kernel`](#mk_spatial_smear_mom_kernel)
   - [`mk_sphere_sum_mom_kernel`](#mk_sphere_sum_mom_kernel)
5. [Fourier-Based Operations](#fourier-based-operations)
   - [`sphere_sum_field`](#sphere_sum_field)
   - [`smear_field`](#smear_field)
   - [`field_convolution`](#field_convolution)
6. [Examples](#examples)

---

## Overview

`field_analysis` provides utilities for analysing and transforming lattice
fields using Fourier-transform techniques. The main capabilities are:

- **Gaussian smearing** — smear a density or observable field in
  coordinate space via multiplication in momentum space. Both 4D and
  spatial-only variants are supported.
- **Sphere summation** — sum field values within a Euclidean-distance
  radius using an FFT-based convolution.
- **Field convolution** — compute the cross-correlation of two complex
  fields via FFT.
- **Step-based smearing** — iterative nearest-neighbour smearing that
  converges to a Gaussian for many steps. Both a local (single-node)
  and a distributed variant are provided.

All Fourier operations use `qlat.mk_fft` and exploit the convolution
theorem: convolution in coordinate space equals pointwise multiplication
in momentum space.

---

## Shift Index Utilities

### `mk_shift_xg_idx_arr`

```python
mk_shift_xg_idx_arr(total_site: tuple, xg_shift: tuple) -> np.ndarray
```

Build an index array that maps each site to the site shifted by
`xg_shift` (with periodic boundary conditions). Only works when the
field lives on a single MPI node (no communication).

| Parameter | Type | Description |
|---|---|---|
| `total_site` | `tuple[int, ...]` | Lattice dimensions, e.g. `(4, 4, 4, 8)` |
| `xg_shift` | `tuple[int, ...]` | Shift in each direction |

Returns an integer `np.ndarray` of shape `(total_volume,)` suitable for
advanced indexing: `field[shift_idx]` yields the shifted field.

Results are cached (up to 16 entries).

---

## Step-Based Smearing

### `smear_field_step_local`

```python
smear_field_step_local(field: Field, coef: float, n_steps: int = 1) -> Field
```

Iterative nearest-neighbour smearing using local index arrays (no MPI
communication). At each step, the field value at each site is mixed with
the average of its 8 nearest neighbours (positive and negative shifts in
all 4 directions).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `field` | `Field` | — | Input field (at least complex type) |
| `coef` | `float` | — | Smearing coefficient (0 = no smearing, 1 = full average) |
| `n_steps` | `int` | `1` | Number of smearing iterations |

Returns a new field; the input is not modified.

**Note:** Requires `geo.num_node() == 1`.

### `smear_field_step`

```python
smear_field_step(field: Field, coef: float, n_steps: int = 1) -> Field
```

Distributed variant of nearest-neighbour smearing. Uses
`field.shift()` which supports MPI communication. The smearing
kernel is identical to `smear_field_step_local`.

Returns a new field; the input is not modified.

---

## Momentum-Space Kernels

### `mk_smear_mom_kernel`

```python
mk_smear_mom_kernel(total_site: tuple, radius: float) -> FieldRealD
```

Build the 4D Gaussian smearing kernel in momentum space:

$$G(k) = \exp\!\Bigl(-\frac{r^2}{4}\sum_{\mu=0}^{3}\sin^2\!\frac{k_\mu}{2}\Bigr)$$

where $k_\mu = 2\pi n_\mu / L_\mu$.

If `radius == np.inf`, the kernel is a delta function at $k = 0$.

Results are cached (up to 128 entries).

### `mk_spatial_smear_mom_kernel`

```python
mk_spatial_smear_mom_kernel(total_site: tuple, radius: float) -> FieldRealD
```

Build the spatial-only Gaussian smearing kernel:

$$G(k) = \exp\!\Bigl(-\frac{2r^2}{3}\sum_{i=0}^{2}\sin^2\!\frac{k_i}{2}\Bigr)$$

The temporal direction is not smeared. If `radius == np.inf`, the
kernel selects the spatial zero mode.

Results are cached (up to 128 entries).

### `mk_sphere_sum_mom_kernel`

```python
mk_sphere_sum_mom_kernel(
    total_site: tuple,
    radius: float,
    is_only_spatial: bool,
) -> FieldRealD
```

Build the momentum-space kernel for sphere summation. In coordinate
space this corresponds to a sharp cutoff: sum all sites within
Euclidean distance `radius`. The kernel is constructed by FFT of a
real-space indicator function.

If `radius == np.inf`, every site contributes equally.

Results are cached (up to 128 entries).

---

## Fourier-Based Operations

### `sphere_sum_field`

```python
sphere_sum_field(
    field: Field,
    radius: float,
    *,
    is_only_spatial: bool = False,
) -> Field
```

Approximate sphere summation of a field via FFT convolution:

$$f_{\text{sphere-summed}}(x) \approx \sum_y \theta(|x-y| < r)\, f(y)$$

The field must be at least complex type. The result is normalised so
that the sum approximates a spherical average.

### `smear_field`

```python
smear_field(
    field: Field,
    radius: float,
    *,
    is_only_spatial: bool = False,
) -> Field
```

Gaussian-smear a field in momentum space and return the result.

$$f_{\text{smear}}(x) \approx
\frac{\sum_y f(y)\,\exp\!\bigl(-(x-y)^2/(2r^2)\bigr)}
     {\sum_y \exp\!\bigl(-y^2/(2r^2)\bigr)}$$

The field must be at least complex type. If `is_only_spatial` is
`True`, only the three spatial directions are smeared.

### `field_convolution`

```python
field_convolution(
    f1: FieldComplexD,
    f2: FieldComplexD,
    idx1: np.ndarray = None,
    idx2: np.ndarray = None,
    *,
    is_only_spatial: bool = False,
) -> FieldComplexD
```

Compute the cross-correlation of two complex fields via FFT:

$$\text{ff}[x_g] = \sum_{x} f_2(x + x_g,\, \text{idx2}) \cdot f_1(x,\, \text{idx1})$$

| Parameter | Type | Default | Description |
|---|---|---|---|
| `f1` | `FieldComplexD` | — | First field |
| `f2` | `FieldComplexD` | — | Second field (must have same `total_site`) |
| `idx1` | `np.ndarray` or `None` | `None` | Multiplicity indices for `f1`; defaults to all |
| `idx2` | `np.ndarray` or `None` | `None` | Multiplicity indices for `f2`; defaults to all |
| `is_only_spatial` | `bool` | `False` | If `True`, convolve only in spatial directions |

`idx1` and `idx2` must have the same length. The result has
multiplicity `len(idx1)`.

---

## Examples

### Gaussian Smearing via FFT

```python
import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

f = q.FieldRealD(geo, 1)
f.set_zero()
f[0] = np.array([1.0])

smeared = q.smear_field(f, radius=2.0)
print(f"Smeared field sum = {np.asarray(smeared).sum():.6f}")

q.end_with_mpi()
```

### Sphere Summation

```python
import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

f = q.FieldRealD(geo, 1)
f.set_zero()
f[0] = np.array([1.0])

summed = q.sphere_sum_field(f, radius=3.0)
print(f"Sphere-summed field sum = {np.asarray(summed).sum():.6f}")

q.end_with_mpi()
```

### Step-Based Smearing (Single Node)

```python
import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

f = q.FieldComplexD(geo, 1)
f.set_zero()
f[0] = np.array([1.0 + 0.0j])

smeared = q.smear_field_step_local(f, coef=0.5, n_steps=10)
print(f"Step-smeared field sum = {np.asarray(smeared).sum():.6f}")

q.end_with_mpi()
```

### Field Convolution

```python
import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

f1 = q.FieldComplexD(geo, 1)
f1.set_zero()
f1[0] = np.array([1.0 + 0.0j])

f2 = q.FieldComplexD(geo, 1)
f2.set_zero()
f2[0] = np.array([1.0 + 0.0j])

ff = q.field_convolution(f1, f2)
print(f"Convolution result at origin = {np.asarray(ff)[0]}")

q.end_with_mpi()
```
