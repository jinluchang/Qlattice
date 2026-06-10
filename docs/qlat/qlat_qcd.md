# `qlat.qcd` — Gauge Fields, Gauge Transformations, and DWF QED Solvers

Source: `qlat/qlat/qcd.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`GaugeField` Class](#gaugefield-class)
   - [Constructor](#constructor)
   - [I/O](#io)
   - [Random Generation](#random-generation)
   - [Unitarization](#unitarization)
   - [Observables](#observables)
   - [Boundary Conditions](#boundary-conditions)
3. [`GaugeTransform` Class](#gaugetransform-class)
   - [Constructor](#constructor-1)
   - [I/O](#io-1)
   - [Random Generation and Unitarization](#random-generation-and-unitarization)
   - [Gauge Transformation Multiplication](#gauge-transformation-multiplication)
   - [Inverse](#inverse)
4. [Gauge Field Functions](#gauge-field-functions)
   - [Plaquette and Link Trace](#plaquette-and-link-trace)
   - [Wilson Lines](#wilson-lines)
   - [Wilson Loop](#wilson-loop)
   - [Random Color Matrix Field](#random-color-matrix-field)
   - [Boundary and Reduction](#boundary-and-reduction)
5. [DWF QED Inverters](#dwf-qed-inverters)
   - [`invert_dwf_qed`](#invert_dwf_qed)
   - [`cg_with_m_dwf_qed`](#cg_with_m_dwf_qed)
   - [`multiply_m_dwf_qed`](#multiply_m_dwf_qed)
6. [Examples](#examples)

---

## Overview

`qcd` provides the core gauge-field infrastructure for lattice QCD
simulations in qlat. It defines two primary classes:

- **`GaugeField`** — A `FieldColorMatrix` with multiplicity 4 (one SU(3)
  matrix per lattice direction) representing the gauge link variables U(x,mu).
  Supports NERSC-format I/O, plaquette and link-trace observables, Wilson
  lines, and boundary condition manipulation.

- **`GaugeTransform`** — A `FieldColorMatrix` with multiplicity 1 representing
  a gauge transformation matrix V(x). Supports application to gauge fields,
  propagators, and fermion fields via the `*` operator, as well as CPS-format
  I/O.

The module also provides standalone functions for computing plaquette fields,
Wilson lines/loops, and domain-wall fermion (DWF) QED inversions that are
used in lattice QED calculations (e.g., for the muon g-2 HLbL contribution).

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf = q.qcd.GaugeField(geo)
rng = q.RngState("test")
gf.set_rand(rng)
print(f"plaquette = {gf.plaq():.10f}")
print(f"link_trace = {gf.link_trace():.10f}")

q.end_with_mpi()
```

---

## `GaugeField` Class

`GaugeField` inherits from `FieldColorMatrix` with fixed multiplicity 4.
Each site stores 4 color matrices — one for each lattice direction mu=0,1,2,3.

### Constructor

#### `GaugeField(geo: Geometry = None, multiplicity: int = 4)`

Create a gauge field. The `multiplicity` argument must be 4 (or 0 for an
uninitialized field).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `geo` | `Geometry` or `None` | `None` | Lattice geometry; `None` creates an empty field |
| `multiplicity` | `int` | `4` | Must be 4 |

---

### I/O

#### `save(path: str)`

Save the gauge field to disk in the standard NERSC format.

#### `load(path: str) -> int`

Load the gauge field from a NERSC-format file. Returns 0 on success.

---

### Random Generation

#### `set_rand(rng: RngState, sigma: float = 0.5, n_step: int = 1)`

Fill the field with random SU(3) matrices generated via Gaussian
smearing around the identity. Parameters `sigma` and `n_step` control the
smearing width and number of smearing steps.

---

### Unitarization

#### `unitarize()`

Project each link matrix to SU(3) in place.

---

### Observables

#### `plaq() -> float`

Return the average plaquette of the gauge field.

#### `link_trace() -> float`

Return the average link trace of the gauge field.

#### `show_info()`

Display the plaquette and link trace to the log.

---

### Boundary Conditions

#### `twist_boundary_at_boundary(lmom: float = -0.5, mu: int = 3)`

Apply a twist boundary condition in direction `mu` by modifying the gauge
field in place. The parameter `lmom` controls the boundary phase (default
-0.5 corresponds to periodic boundary conditions).

---

## `GaugeTransform` Class

`GaugeTransform` inherits from `FieldColorMatrix` with fixed multiplicity 1.
Each site stores one SU(3) matrix V(x) representing a gauge transformation.

### Constructor

#### `GaugeTransform(geo: Geometry = None, multiplicity: int = 1)`

Create a gauge transformation field. The `multiplicity` argument must be 1
(or 0 for an uninitialized field).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `geo` | `Geometry` or `None` | `None` | Lattice geometry |
| `multiplicity` | `int` | `1` | Must be 1 |

---

### I/O

#### `save(path: str)`

Save as RealD precision using the generic field format (`save_double`).

#### `load(path: str) -> int`

Load as RealD precision using the generic field format (`load_double`).

#### `save_cps(path: str)`

Save in the CPS format.

#### `load_cps(path: str) -> int`

Load from the CPS format.

---

### Random Generation and Unitarization

#### `set_rand(rng: RngState, sigma: float = 0.5, n_step: int = 1)`

Fill with random SU(3) matrices via Gaussian smearing.

#### `unitarize()`

Project each matrix to SU(3) in place.

---

### Gauge Transformation Multiplication

#### `__mul__(other) -> GaugeTransform | GaugeField | Prop | SelProp | PselProp | FermionField4d | list`

Apply the gauge transformation to another object. The result type depends
on the operand:

| Operand type | Result type | Operation |
|---|---|---|
| `GaugeTransform` | `GaugeTransform` | Compose transformations: V * V' |
| `GaugeField` | `GaugeField` | Transform gauge links: V * U |
| `Prop` | `Prop` | Transform full propagator |
| `SelProp` | `SelProp` | Transform selected-point propagator |
| `PselProp` | `PselProp` | Transform point-selected propagator |
| `FermionField4d` | `FermionField4d` | Transform fermion field |
| `list` | `list` | Element-wise transformation |

---

### Inverse

#### `inv() -> GaugeTransform`

Return the inverse transformation V^dagger.

---

## Gauge Field Functions

### Plaquette and Link Trace

#### `gf_avg_plaq(gf: GaugeField) -> float`

Return the average plaquette.

#### `gf_avg_spatial_plaq(gf: GaugeField) -> float`

Return the average spatial plaquette (mu,nu with mu,nu in {0,1,2}).

#### `gf_avg_link_trace(gf: GaugeField) -> float`

Return the average link trace (trace of U divided by Nc).

#### `gf_show_info(gf: GaugeField)`

Display the plaquette and link trace to the log.

#### `gf_plaq_field(gf: GaugeField) -> FieldRealD`

Return a `FieldRealD` with multiplicity 1 containing the plaquette value
at each site.

#### `set_tr_less_anti_herm_matrix(fc: FieldColorMatrix)`

Make each matrix in the field traceless anti-Hermitian in place.

---

### Wilson Lines

#### `gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n=None)`

Compute a Wilson line along `path` and store the result in multiplicity `m`
of `wlf`. The gauge field `gf_ext` must already include halo expansion.

| Parameter | Type | Description |
|---|---|---|
| `wlf` | `FieldColorMatrix` | Output field (modified in multiplicity `m`) |
| `m` | `int` | Multiplicity index to write |
| `gf_ext` | `GaugeField` | Expanded gauge field |
| `path` | `list[int]` | Sequence of directions; negative values mean reversed links |
| `path_n` | `list[int]` or `None` | Optional repetition counts per step |

Path encoding: `mu` for forward link, `-mu-1` for backward link.
For example, `[0, 0, 3, -1, -1]` means mu, mu, nu, -mu, -mu.

#### `gf_wilson_lines_no_comm(gf_ext, path_list) -> FieldColorMatrix`

Compute multiple Wilson lines. Each entry in `path_list` can be a plain
`path` list or a `(path, path_n)` tuple.

| Parameter | Type | Description |
|---|---|---|
| `gf_ext` | `GaugeField` | Expanded gauge field |
| `path_list` | `list` | List of path specifications |

Returns a `FieldColorMatrix` with multiplicity `len(path_list)`.

---

### Wilson Loop

#### `gf_avg_wilson_loop_normalized_tr(gf: GaugeField, l: int, t: int) -> float`

Return the average Wilson loop trace (normalized by Nc=3) for an
l x t rectangular loop.

| Parameter | Type | Description |
|---|---|---|
| `gf` | `GaugeField` | Gauge field |
| `l` | `int` | Spatial extent of the loop |
| `t` | `int` | Temporal extent of the loop |

---

### Random Color Matrix Field

#### `set_g_rand_color_matrix_field(fc: FieldColorMatrix, rng: RngState, sigma: float, n_steps: int = 1)`

Fill `fc` with random SU(3) matrices generated via Gaussian smearing
around the identity.

---

### Boundary and Reduction

#### `gf_twist_boundary_at_boundary(gf: GaugeField, lmom: float = -0.5, mu: int = 3)`

Apply a twist boundary condition in direction `mu` by modifying `gf` in place.

#### `mk_left_expanded_field(gf) -> FieldColorMatrix`

Return a left-expanded copy of the gauge field (expansion of 1 in each
direction on the left side). Equivalent to `set_left_expanded_gauge_field`
in the C++ library.

#### `gf_reduce_half(gf: GaugeField) -> GaugeField`

Return a new gauge field with half the lattice size in all directions.
Does not modify the input. Can be combined with shifts, e.g.,
`gf_reduce_half(gf.shift())`.

---

## DWF QED Inverters

These functions solve domain-wall fermion (DWF) linear systems with
QED corrections. They require an expanded gauge field `gf1` (obtained via
`q.mk_left_expanded_field(gf)`).

The optional `t_wick_phase_factor_arr` parameter provides per-timeslice
Wick rotation phase factors for QED photon propagators.

### `invert_dwf_qed`

```python
invert_dwf_qed(f_in4d, gf1, mass, m5, ls, *, t_wick_phase_factor_arr=None,
               is_dagger=False, stop_rsd=1e-8, max_num_iter=50000) -> FieldComplexD
```

Solve `M * out = in` (or `M^dag * out = in` if `is_dagger=True`) where M
is the DWF operator with QED corrections. The input `f_in4d` is a 4-D
fermion field and the output is also 4-D (properly projected from the 5-D
solution).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `f_in4d` | `FieldComplexD` | — | 4-D source field |
| `gf1` | `FieldComplexD` | — | Left-expanded gauge field |
| `mass` | `float` | — | Quark mass |
| `m5` | `float` | — | Domain-wall height |
| `ls` | `int` | — | Fifth-dimension extent |
| `t_wick_phase_factor_arr` | array or `None` | `None` | Per-timeslice Wick phase factors |
| `is_dagger` | `bool` | `False` | Use M^dag if True |
| `stop_rsd` | `float` | `1e-8` | Stopping residual |
| `max_num_iter` | `int` | `50000` | Maximum CG iterations |

### `cg_with_m_dwf_qed`

```python
cg_with_m_dwf_qed(f_in5d, gf1, mass, m5, ls, *, t_wick_phase_factor_arr=None,
                  is_dagger=False, stop_rsd=1e-8, max_num_iter=50000) -> FieldComplexD
```

Solve `M^dag M * out = in` (or `M M^dag * out = in` if `is_dagger=True`).
Input and output are 5-D fermion fields.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `f_in5d` | `FieldComplexD` | — | 5-D source field |
| `gf1` | `FieldComplexD` | — | Left-expanded gauge field |
| `mass` | `float` | — | Quark mass |
| `m5` | `float` | — | Domain-wall height |
| `ls` | `int` | — | Fifth-dimension extent |
| `t_wick_phase_factor_arr` | array or `None` | `None` | Per-timeslice Wick phase factors |
| `is_dagger` | `bool` | `False` | Use M M^dag if True |
| `stop_rsd` | `float` | `1e-8` | Stopping residual |
| `max_num_iter` | `int` | `50000` | Maximum CG iterations |

### `multiply_m_dwf_qed`

```python
multiply_m_dwf_qed(f_in5d, gf1, mass, m5, ls, *, t_wick_phase_factor_arr=None,
                   is_dagger=False) -> FieldComplexD
```

Apply the DWF operator with QED corrections: compute `out = M * in` (or
`out = M^dag * in` if `is_dagger=True`). Input and output are 5-D fields.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `f_in5d` | `FieldComplexD` | — | 5-D input field |
| `gf1` | `FieldComplexD` | — | Left-expanded gauge field |
| `mass` | `float` | — | Quark mass |
| `m5` | `float` | — | Domain-wall height |
| `ls` | `int` | — | Fifth-dimension extent |
| `t_wick_phase_factor_arr` | array or `None` | `None` | Per-timeslice Wick phase factors |
| `is_dagger` | `bool` | `False` | Apply M^dag if True |

---

## Examples

### Basic Gauge Field Operations

```python
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf = q.qcd.GaugeField(geo)

rng = q.RngState("test-seed")
gf.set_rand(rng, sigma=0.5, n_step=3)
gf.unitarize()

print(f"plaq      = {gf.plaq():.10f}")
print(f"link_tr   = {gf.link_trace():.10f}")
gf.show_info()

q.end_with_mpi()
```

### Save and Load Gauge Field (NERSC Format)

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf = q.qcd.GaugeField(geo)

rng = q.RngState("save-load-test")
gf.set_rand(rng)

gf.save("gauge_field.nersc")
gf2 = q.qcd.GaugeField()
gf2.load("gauge_field.nersc")

assert abs(gf.plaq() - gf2.plaq()) < 1e-12

q.end_with_mpi()
```

### Gauge Transformation

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

gf = q.qcd.GaugeField(geo)
rng = q.RngState("gt-test")
gf.set_rand(rng, sigma=0.3, n_step=3)
gf.unitarize()

plaq_before = gf.plaq()

gt = q.qcd.GaugeTransform(geo)
gt.set_rand(rng, sigma=0.1, n_step=1)
gt.unitarize()

gf_transformed = gt * gf
print(f"plaq before transform = {plaq_before:.10f}")
print(f"plaq after  transform = {gf_transformed.plaq():.10f}")

gt_inv = gt.inv()
gf_roundtrip = gt_inv * gf_transformed
print(f"plaq after inverse    = {gf_roundtrip.plaq():.10f}")

q.end_with_mpi()
```

### Wilson Lines

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

gf = q.qcd.GaugeField(geo)
rng = q.RngState("wilson-line")
gf.set_rand(rng, sigma=0.3, n_step=3)
gf.unitarize()

gf_ext = q.mk_left_expanded_field(gf)

mu, nu = 0, 1
path_list = [
    [mu, mu, nu, -mu - 1, -mu - 1],
    ([mu, nu, -mu - 1], [2, 1, 2]),
]
wlf = q.qcd.gf_wilson_lines_no_comm(gf_ext, path_list)
print(f"Wilson line field multiplicity = {wlf.multiplicity}")

q.end_with_mpi()
```

### DWF QED Inversion

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

gf = q.qcd.GaugeField(geo)
rng = q.RngState("dwf-qed")
gf.set_rand(rng, sigma=0.3, n_step=3)
gf.unitarize()

gf1 = q.mk_left_expanded_field(gf)

mass = 0.01
m5 = 1.8
ls = 12

f_in = q.FieldComplexD(geo, 1)
f_in.set_rand(q.RngState("source"))

f_out = q.qcd.invert_dwf_qed(f_in, gf1, mass, m5, ls, stop_rsd=1e-8)
print(f"DWF QED inversion completed")

q.end_with_mpi()
```
