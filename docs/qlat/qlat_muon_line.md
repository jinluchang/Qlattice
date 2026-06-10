# `qlat.muon_line` — Muon-Line Integrals for HLbL Scattering

Source: `qlat/qlat/muon_line.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [CUBA Integration Library](#cuba-integration-library)
3. [Interpolation Management](#interpolation-management)
   - [`clear_muon_line_interpolations`](#clear_muon_line_interpolations)
   - [`get_number_of_muon_line_interpolations`](#get_number_of_muon_line_interpolations)
   - [`compute_save_muonline_interpolation`](#compute_save_muonline_interpolation)
   - [`load_multiple_muonline_interpolations`](#load_multiple_muonline_interpolations)
4. [Extra Weights](#extra-weights)
   - [`get_muon_line_m_extra_weights`](#get_muon_line_m_extra_weights)
   - [`set_muon_line_m_extra_weights`](#set_muon_line_m_extra_weights)
5. [Core Computation Functions](#core-computation-functions)
   - [`calc_muon_line_m`](#calc_muon_line_m)
   - [`get_muon_line_m`](#get_muon_line_m)
   - [`get_muon_line_m_extra`](#get_muon_line_m_extra)
   - [`get_muon_line_m_extra_lat`](#get_muon_line_m_extra_lat)
6. [Examples](#examples)

---

## Overview

`muon_line` computes the muon-line tensor

$$
\mathcal{M}_{i,\rho,\sigma,\lambda}(x,y,z)
$$

that appears in the hadronic light-by-light (HLbL) scattering contribution
to the muon anomalous magnetic moment :math:`(g-2)_\mu`.  The computation
follows the method described in [arXiv:2304.04423](https://arxiv.org/abs/2304.04423), Eq.~(9).

The module provides:

- **Direct numerical integration** via `calc_muon_line_m` (uses the CUBA
  multi-dimensional integration library).
- **Pre-computed interpolation tables** that can be saved to disk and loaded
  later, avoiding repeated expensive integrals.
- **Extra-weighted variants** (`get_muon_line_m_extra`,
  `get_muon_line_m_extra_lat`) that combine multiple interpolation tables
  with configurable weights for subtracted and unsubtracted contributions.

All core functions return a 4-D NumPy array of shape `(3, 4, 4, 4)` and
dtype `float64`, where the indices are `[i, rho, sigma, lambda]`.

The module sets the environment variable `CUBACORES=0` on import to
control CUBA's threading behaviour.

---

## CUBA Integration Library

### `has_cuba() -> bool`

Return `True` if the CUBA integration library is available at runtime.

```python
import qlat as q

q.begin_with_mpi(size_node_list)

if q.has_cuba():
    print("CUBA is available")

q.end_with_mpi()
```

### `test_integration_multi_dimensional()`

Run a built-in test of the CUBA multi-dimensional integration routines.
Useful for verifying that the integration library works correctly.

---

## Interpolation Management

Interpolations are pre-computed tables of the muon-line tensor stored on
disk.  Once loaded they allow fast retrieval of
:math:`\mathcal{M}_{i,\rho,\sigma,\lambda}` without re-evaluating the
integral.

### `clear_muon_line_interpolations()`

Release all loaded interpolation tables from memory.

### `get_number_of_muon_line_interpolations() -> int`

Return the number of interpolation tables currently loaded in memory.

### `compute_save_muonline_interpolation(path, dims, eps) -> int`

Compute and save a muon-line interpolation table to disk.

| Parameter | Type | Description |
|---|---|---|
| `path` | `str` | Directory path where the interpolation is saved |
| `dims` | `list[int]` | Grid resolution in each of the 5 dimensions, e.g. `[6, 6, 6, 6, 6]` up to `[16, 16, 16, 16, 16]` |
| `eps` | `tuple` | Integration parameters: `(epsabs, epsrel, mineval, maxeval)` |

Default integration tolerances: `epsabs=1e-8`, `epsrel=1e-3`.

### `load_multiple_muonline_interpolations(path, idx_list) -> int`

Load pre-computed interpolation tables from `path`.

| Parameter | Type | Description |
|---|---|---|
| `path` | `str` | Root directory containing interpolation checkpoints |
| `idx_list` | `list[int]` | Indices to load.  If empty, load all available tables |

Tables are expected under `path/{idx:010d}/checkpoint`.  Returns the number
of tables loaded.

---

## Extra Weights

Extra weights control how multiple interpolation tables are combined in
`get_muon_line_m_extra` and `get_muon_line_m_extra_lat`.  The weights are
a list of lists of floats, one inner list per interpolation index.

### `get_muon_line_m_extra_weights() -> list[list[float]]`

Return the current extra weights.

### `set_muon_line_m_extra_weights(weights=None)`

Set the extra weights.  Pass `None` to restore the default values.

---

## Core Computation Functions

All functions below return a NumPy array of shape `(3, 4, 4, 4)` with
dtype `float64`.

### `calc_muon_line_m`

```python
calc_muon_line_m(x: CoordinateD, y: CoordinateD, eps: tuple) -> numpy.ndarray
```

Directly compute the muon-line tensor by numerical integration (CUBA).

| Parameter | Type | Description |
|---|---|---|
| `x` | `CoordinateD` | First spacetime coordinate |
| `y` | `CoordinateD` | Second spacetime coordinate |
| `eps` | `tuple` | `(epsabs, epsrel, mineval, maxeval)` |

Wraps the C++ function `muon_line_sym_py`.

### `get_muon_line_m`

```python
get_muon_line_m(x: CoordinateD, y: CoordinateD, z: CoordinateD,
                idx: int, eps: tuple) -> numpy.ndarray
```

Return the muon-line tensor
:math:`\mathcal{M}_{i,\rho,\sigma,\lambda}(x,y,z)`.

| Parameter | Type | Description |
|---|---|---|
| `x` | `CoordinateD` | First spacetime coordinate |
| `y` | `CoordinateD` | Second spacetime coordinate |
| `z` | `CoordinateD` | Third spacetime coordinate |
| `idx` | `int` | Interpolation index.  If `idx < 0`, compute directly via integration (ignores loaded tables).  Otherwise use the loaded interpolation at `idx`. |
| `eps` | `tuple` | `(epsabs, epsrel, mineval, maxeval)` — only used when `idx < 0` |

Wraps the C++ function `get_muon_line_m_py`.

### `get_muon_line_m_extra`

```python
get_muon_line_m_extra(x: CoordinateD, y: CoordinateD, z: CoordinateD,
                      tag: int) -> numpy.ndarray
```

Return the muon-line tensor using weighted combinations of loaded
interpolations.

| Parameter | Type | Description |
|---|---|---|
| `x` | `CoordinateD` | First spacetime coordinate |
| `y` | `CoordinateD` | Second spacetime coordinate |
| `z` | `CoordinateD` | Third spacetime coordinate |
| `tag` | `int` | `0` for subtracted, `1` for unsubtracted |

The following interpolation tables must be loaded:

| Index | Grid | Type | Used by `tag` |
|---|---|---|---|
| 0 | 6^5 | with-sub | |
| 1 | 8^5 | with-sub | `tag=0` |
| 2 | 10^5 | with-sub | |
| 3 | 12^5 | with-sub | `tag=0` |
| 4 | 14^5 | with-sub | |
| 5 | 16^5 | with-sub | `tag=0` |
| 6 | 6^5 | no-sub | |
| 7 | 8^5 | no-sub | `tag=1` |
| 8 | 10^5 | no-sub | |
| 9 | 12^5 | no-sub | `tag=1` |
| 10 | 14^5 | no-sub | |
| 11 | 16^5 | no-sub | `tag=1` |

Wraps the C++ function `get_muon_line_m_extra_py`.

### `get_muon_line_m_extra_lat`

```python
get_muon_line_m_extra_lat(x: Coordinate, y: Coordinate, z: Coordinate,
                          total_site: Coordinate, a: float,
                          tag: int) -> numpy.ndarray
```

Lattice-aware version of `get_muon_line_m_extra`.  Coordinates are given
as integer lattice sites and `a` is the muon mass in lattice units.

| Parameter | Type | Description |
|---|---|---|
| `x` | `Coordinate` | First lattice site |
| `y` | `Coordinate` | Second lattice site |
| `z` | `Coordinate` | Third lattice site |
| `total_site` | `Coordinate` | Total lattice dimensions |
| `a` | `float` | Muon mass in lattice units |
| `tag` | `int` | `0` for subtracted, `1` for unsubtracted |

Requires the same set of interpolation tables as `get_muon_line_m_extra`.

Wraps the C++ function `get_muon_line_m_extra_lat_py`.

---

## Examples

### Check CUBA Availability

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

print("CUBA available:", q.has_cuba())

q.end_with_mpi()
```

### Compute Muon-Line Tensor Directly

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

x = q.CoordinateD([0.1, 0.2, 0.3, 0.4])
y = q.CoordinateD([0.5, 0.6, 0.7, 0.8])

eps = (1e-8, 1e-3, 0, 0)
m = q.calc_muon_line_m(x, y, eps)

print("Shape:", m.shape)       # (3, 4, 4, 4)
print("dtype:", m.dtype)       # float64

q.end_with_mpi()
```

### Load Interpolations and Query

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

# Load pre-computed interpolation tables
path = "/path/to/interpolations"
q.load_multiple_muonline_interpolations(path, [])
print("Loaded:", q.get_number_of_muon_line_interpolations())

# Use a loaded interpolation (idx >= 0)
x = q.CoordinateD([0.1, 0.2, 0.3, 0.4])
y = q.CoordinateD([0.5, 0.6, 0.7, 0.8])
z = q.CoordinateD([0.0, 0.0, 0.0, 0.0])
m = q.get_muon_line_m(x, y, z, idx=1, eps=(0, 0, 0, 0))

print("Shape:", m.shape)

# Clean up
q.clear_muon_line_interpolations()

q.end_with_mpi()
```

### Use Extra-Weighted Muon-Line Tensor (Lattice Version)

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

# Load all required interpolation tables (indices 0-11)
path = "/path/to/interpolations"
q.load_multiple_muonline_interpolations(path, [])

# Integer lattice coordinates
x = q.Coordinate([1, 0, 0, 0])
y = q.Coordinate([0, 1, 0, 0])
z = q.Coordinate([0, 0, 0, 0])
total_site = q.Coordinate([16, 16, 16, 32])
a = 0.01  # muon mass in lattice units

# Subtracted contribution (tag=0)
m_sub = q.get_muon_line_m_extra_lat(x, y, z, total_site, a, tag=0)

# Unsubtracted contribution (tag=1)
m_nosub = q.get_muon_line_m_extra_lat(x, y, z, total_site, a, tag=1)

print("sub shape:", m_sub.shape)      # (3, 4, 4, 4)
print("nosub shape:", m_nosub.shape)  # (3, 4, 4, 4)

q.clear_muon_line_interpolations()

q.end_with_mpi()
```
