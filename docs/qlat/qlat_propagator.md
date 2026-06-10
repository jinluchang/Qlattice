# `qlat.propagator` — Propagator Types and Utilities

Source: `qlat/qlat/propagator.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Propagator Classes](#propagator-classes)
   - [`Prop`](#prop)
   - [`SelProp`](#selprop)
   - [`PselProp`](#pselprop)
   - [`SpinProp`](#spinprop)
   - [`FermionField4d`](#fermionfield4d)
3. [Source Construction](#source-construction)
   - [`set_point_src`](#set_point_src)
   - [`set_wall_src`](#set_wall_src)
   - [`mk_point_src`](#mk_point_src)
   - [`mk_wall_src`](#mk_wall_src)
   - [`set_rand_vol_u1` / `mk_rand_vol_u1`](#set_rand_vol_u1--mk_rand_vol_u1)
   - [`set_rand_vol_u1_src` / `mk_rand_vol_u1_src`](#set_rand_vol_u1_src--mk_rand_vol_u1_src)
   - [`mk_rand_u1_src`](#mk_rand_u1_src)
   - [`get_rand_u1_sol`](#get_rand_u1_sol)
   - [`mk_rand_u1_prop`](#mk_rand_u1_prop)
4. [Inversion Functions](#inversion-functions)
   - [`free_invert`](#free_invert)
   - [`invert_qed`](#invert_qed)
5. [Format Conversion](#format-conversion)
   - [`convert_mspincolor_from_wm`](#convert_mspincolor_from_wm)
   - [`convert_wm_from_mspincolor`](#convert_wm_from_mspincolor)
6. [Fermion Field Utilities](#fermion-field-utilities)
   - [`mk_ff_list_from_prop`](#mk_ff_list_from_prop)
   - [`mk_prop_from_ff_list`](#mk_prop_from_ff_list)
7. [Miscellaneous](#miscellaneous)
   - [`flip_tpbc_with_tslice`](#flip_tpbc_with_tslice)
   - [`free_scalar_invert_mom_cfield`](#free_scalar_invert_mom_cfield)
   - [`free_scalar_invert_cfield`](#free_scalar_invert_cfield)
8. [Examples](#examples)

---

## Overview

`propagator` defines the primary propagator container types used in
domain-wall fermion lattice QCD calculations and provides functions for
source construction, inversion, and format conversion.

The five propagator/field types are:

| Type | Element Type | Description |
|---|---|---|
| `Prop` | `WilsonMatrix` (12x12) | Full-volume propagator field |
| `SelProp` | `WilsonMatrix` | Propagator on a `FieldSelection` subset of sites |
| `PselProp` | `WilsonMatrix` | Propagator on a `PointsSelection` subset of sites |
| `SpinProp` | `SpinMatrix` (4x4) | Full-volume spin propagator (no colour index) |
| `FermionField4d` | `WilsonVector` (12-component) | Single fermion field |

Each class enforces `multiplicity == 1` — there is exactly one matrix or
vector element per site.

`Prop` stores a 12x12 complex matrix (4 spin x 3 colour) at every lattice
site.  `SelProp` and `PselProp` are selected-site variants used for
point-source solves and measurements.  `SpinProp` is the spin-only (4x4)
counterpart used in QED corrections.  `FermionField4d` is a single
12-component Wilson vector field.

---

## Propagator Classes

### `Prop`

```python
class Prop(FieldWilsonMatrix)
```

Full-volume propagator.  Each site stores one `WilsonMatrix` (12x12
complex matrix, multiplicity = 1).

#### Constructors

| Signature | Description |
|---|---|
| `Prop()` | Uninitialized propagator |
| `Prop(geo: Geometry)` | Allocate on the given geometry |

#### Methods

| Method | Description |
|---|---|
| `get_elem_wm(index, m=0) -> WilsonMatrix` | Return the Wilson matrix at site `index` |
| `glb_sum_tslice(t_dir=3) -> PselProp` | Sum propagator over each time slice, returning a `PselProp` |
| `__getstate__` / `__setstate__` | Pickle support (single-node only) |

### `SelProp`

```python
class SelProp(SelectedFieldWilsonMatrix)
```

Propagator defined on a field selection (a subset of lattice sites).

#### Constructors

| Signature | Description |
|---|---|
| `SelProp()` | Uninitialized |
| `SelProp(fsel: FieldSelection)` | Allocate on the selected sites |

#### Methods

| Method | Description |
|---|---|
| `get_elem_wm(idx, m=0) -> WilsonMatrix` | Return the Wilson matrix at selected-site index `idx` |
| `__getstate__` / `__setstate__` | Pickle support (single-node only) |

### `PselProp`

```python
class PselProp(SelectedPointsWilsonMatrix)
```

Propagator defined on a points selection (an ordered list of sites).

#### Constructors

| Signature | Description |
|---|---|
| `PselProp()` | Uninitialized |
| `PselProp(psel: PointsSelection)` | Allocate for the given point set |
| `PselProp(psel, multiplicity)` | Must have `multiplicity == 1` |
| `PselProp(sel_prop: SelProp)` | Copy selected sites from a `SelProp` |

#### Methods

| Method | Description |
|---|---|
| `get_elem_wm(idx, m=0) -> WilsonMatrix` | Return the Wilson matrix at point index `idx` |
| `__getstate__` / `__setstate__` | Pickle support (single-node only) |

### `SpinProp`

```python
class SpinProp(FieldSpinMatrix)
```

Full-volume spin propagator (4x4 complex matrix, no colour index).
Used for QED corrections.

#### Constructors

| Signature | Description |
|---|---|
| `SpinProp()` | Uninitialized |
| `SpinProp(geo: Geometry)` | Allocate on the given geometry |

#### Methods

| Method | Description |
|---|---|
| `get_elem_sm(index, m=0) -> SpinMatrix` | Return the spin matrix at site `index` |
| `glb_sum_tslice(t_dir=3) -> SelectedPointsSpinMatrix` | Sum over each time slice |
| `__getstate__` / `__setstate__` | Pickle support (single-node only) |

### `FermionField4d`

```python
class FermionField4d(FieldWilsonVector)
```

A single 12-component fermion field (one `WilsonVector` per site).

#### Constructors

| Signature | Description |
|---|---|
| `FermionField4d()` | Uninitialized |
| `FermionField4d(geo: Geometry)` | Allocate on the given geometry |

---

## Source Construction

### `set_point_src`

```python
set_point_src(prop_src: Prop, geo: Geometry, xg: Coordinate,
              value: complex = 1.0) -> None
```

Set `prop_src` to a point source at global coordinate `xg` with the given
colour-spin amplitude.

### `set_wall_src`

```python
set_wall_src(prop_src: Prop, geo: Geometry, tslice: int,
             lmom: CoordinateD = None) -> None
```

Set `prop_src` to a wall source at time slice `tslice` with optional
spatial momentum `lmom`.  If `lmom` is `None`, zero momentum is used.

### `mk_point_src`

```python
mk_point_src(geo: Geometry, xg: Coordinate,
             value: complex = 1.0) -> Prop
```

Create and return a new `Prop` with a point source.  Convenience wrapper
around `set_point_src`.

### `mk_wall_src`

```python
mk_wall_src(geo: Geometry, tslice: int,
            lmom: CoordinateD = None) -> Prop
```

Create and return a new `Prop` with a wall source.  Convenience wrapper
around `set_wall_src`.

### `set_rand_vol_u1` / `mk_rand_vol_u1`

```python
set_rand_vol_u1(fu1: FieldComplexD, geo: Geometry,
                multiplicity: int, rs: RngState) -> None

mk_rand_vol_u1(geo: Geometry, multiplicity: int,
               rs: RngState) -> FieldComplexD
```

Generate a volume-filling random U(1) field.  `mk_rand_vol_u1` creates
and returns the field; `set_rand_vol_u1` fills an existing one.

### `set_rand_vol_u1_src` / `mk_rand_vol_u1_src`

```python
set_rand_vol_u1_src(prop_src: Prop, fu1: FieldComplexD) -> None

mk_rand_vol_u1_src(fu1: FieldComplexD) -> Prop
```

Construct a random-U1 volume source from a U(1) field `fu1`.  The source
propagator is proportional to `fu1`.

### `mk_rand_u1_src`

```python
mk_rand_u1_src(sel, rs: RngState) -> tuple[Prop, FieldComplexD]
```

Create a random-U1 source on selected sites.  `sel` can be a
`FieldSelection` or a `PointsSelection`.

Returns `(prop_src, fu1)` where `fu1` stores the random U(1) numbers.

### `get_rand_u1_sol`

```python
get_rand_u1_sol(prop_sol: Prop, fu1: FieldComplexD,
                sel) -> SelProp | PselProp
```

Extract the solution on selected sites after solving with a random-U1
source.  `sel` can be a `FieldSelection` (returns `SelProp`) or a
`PointsSelection` (returns `PselProp`).

### `mk_rand_u1_prop`

```python
mk_rand_u1_prop(inv, sel, rs: RngState) -> SelProp | PselProp
```

High-level interface: generate a random-U1 source, invert, and return the
solution on selected sites.  `inv` is an inverter object (supports
`inv * prop_src`).  `sel` can be a `FieldSelection` or `PointsSelection`.

---

## Inversion Functions

### `free_invert`

```python
free_invert(prop_src, mass: float, m5: float = 1.0,
            momtwist: CoordinateD = None) -> Prop | SpinProp
```

Compute the free (gauge-field-independent) inverse of a propagator source.
Supports both `Prop` and `SpinProp` input; returns the same type.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `prop_src` | `Prop` or `SpinProp` | — | Source field |
| `mass` | `float` | — | Fermion mass |
| `m5` | `float` | `1.0` | Fifth-dimensional mass for DWF |
| `momtwist` | `CoordinateD` | `[0,0,0,0]` | Twist angles for twisted-boundary conditions |

### `invert_qed`

```python
invert_qed(sp_src: SpinProp, gf1: FieldComplexD,
           mass: float, m5: float, ls: int, *,
           t_wick_phase_factor_arr=None,
           is_dagger: bool = False,
           stop_rsd: float = 1e-8,
           max_num_iter: int = 50000) -> SpinProp
```

Invert the Dirac operator in the presence of a QED gauge field for
spin-propagator sources.  Used for electromagnetic corrections.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sp_src` | `SpinProp` | — | Spin-propagator source |
| `gf1` | `FieldComplexD` | — | Left-expanded QED gauge field |
| `mass` | `float` | — | Fermion mass |
| `m5` | `float` | — | DWF fifth-dimensional mass |
| `ls` | `int` | — | Fifth-dimensional extent |
| `t_wick_phase_factor_arr` | array-like | `None` | Optional Wick rotation phase factors per time slice |
| `is_dagger` | `bool` | `False` | Whether to use the conjugate operator |
| `stop_rsd` | `float` | `1e-8` | Stopping residual |
| `max_num_iter` | `int` | `50000` | Maximum CG iterations |

---

## Format Conversion

### `convert_mspincolor_from_wm`

```python
convert_mspincolor_from_wm(prop_wm) -> Prop | SelProp | PselProp
```

Convert a propagator from Wilson-matrix layout (spin-colour interleaved as
a 12x12 matrix) to m-spin-colour layout.  Accepts `Prop`, `SelProp`, or
`PselProp`; returns the same type.

### `convert_wm_from_mspincolor`

```python
convert_wm_from_mspincolor(prop_msc) -> Prop | SelProp | PselProp
```

Convert from m-spin-colour layout back to Wilson-matrix layout.  Inverse
of `convert_mspincolor_from_wm`.

---

## Fermion Field Utilities

### `mk_ff_list_from_prop`

```python
mk_ff_list_from_prop(prop: Prop) -> list[FermionField4d]
```

Split a `Prop` into a list of 12 `FermionField4d` objects (one per
colour-spin column of the Wilson matrix).

### `mk_prop_from_ff_list`

```python
mk_prop_from_ff_list(ff_list: list[FermionField4d]) -> Prop
```

Reassemble a `Prop` from a list of 12 `FermionField4d` objects.  Inverse
of `mk_ff_list_from_prop`.

---

## Miscellaneous

### `flip_tpbc_with_tslice`

```python
flip_tpbc_with_tslice(prop, tslice_flip_tpbc) -> None
```

Flip the temporal boundary condition at `tslice_flip_tpbc` for a `SelProp`
or `PselProp`.

### `free_scalar_invert_mom_cfield`

```python
free_scalar_invert_mom_cfield(f: FieldComplexD, mass: float) -> None
```

Apply the free-scalar inverse in momentum space to a complex field `f`
(modified in place).

### `free_scalar_invert_cfield`

```python
free_scalar_invert_cfield(src: FieldComplexD, mass: float,
                          *, mode_fft: int = 1) -> FieldComplexD
```

Compute the free-scalar inverse of a complex field in position space.
Transforms to momentum space, applies the inverse, and transforms back.

---

## Examples

### Point Source and Free Inversion

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
xg = q.Coordinate([0, 0, 0, 0])

prop_src = q.mk_point_src(geo, xg)
prop_sol = q.free_invert(prop_src, mass=0.1)

wm = prop_sol.get_elem_wm(0)
print("Wilson matrix at site 0 shape:", q.asarray(wm).shape)

q.end_with_mpi()
```

### Wall Source with Momentum

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

lmom = q.CoordinateD([0.0, 0.0, 0.0, 0.0])
prop_src = q.mk_wall_src(geo, tslice=0, lmom=lmom)

print("Wall source created on geometry:", geo)

q.end_with_mpi()
```

### Random U1 Source with Selection

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

# Create a points selection
n_points = 8
psel = q.PointsSelection(geo, [q.Coordinate([0, 0, 0, i]) for i in range(n_points)])

rs = q.RngState("test-seed")
prop_src, fu1 = q.mk_rand_u1_src(psel, rs)

print("Source prop type:", type(prop_src).__name__)
print("U1 field type:", type(fu1).__name__)

q.end_with_mpi()
```

### Splitting and Reassembling a Propagator

```python
import qlat as q

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

prop_src = q.mk_point_src(geo, q.Coordinate([0, 0, 0, 0]))
prop_sol = q.free_invert(prop_src, mass=0.1)

# Split into 12 fermion fields
ff_list = q.mk_ff_list_from_prop(prop_sol)
print(f"Number of fermion fields: {len(ff_list)}")

# Reassemble
prop_rebuilt = q.mk_prop_from_ff_list(ff_list)

q.end_with_mpi()
```

### Pickle Round-Trip

```python
import qlat as q
import pickle

size_node_list = [[1, 1, 1, 1]]

q.begin_with_mpi(size_node_list)

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)

prop = q.mk_point_src(geo, q.Coordinate([0, 0, 0, 0]))

data = pickle.dumps(prop)
prop2 = pickle.loads(data)

print("Pickle round-trip succeeded")

q.end_with_mpi()
```
