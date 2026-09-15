# `qlat.vector_utils` — Legacy Vector Utilities

Source: `qlat/qlat/vector_utils.pyx`

> **Note:** Update this document when updating the source file.

## Outline

- `diff_gauge(g0, g1)` — difference between two gauge fields.
- `diff_prop(p0, p1)` — difference between two propagators.
- `load_gwu_link(g0, path)` — read a gauge field in gwu link format.
- `save_gwu_prop(prop, path)` — write a propagator in gwu format.
- `load_gwu_prop(prop, path)` — read a propagator in gwu format.
- `save_gwu_noiP(prop, path)` — write the `noiP` projection of a propagator.
- `load_gwu_noiP(prop, path)` — read the `noiP` projection of a propagator.
- `load_qlat_link(g0, path)` — read a gauge field in qlat format.
- `save_qlat_link(g0, path)` — write a gauge field in qlat format.
- `random_point_src(prop, seed)` — random point source.
- `make_point_prop(prop, sp)` — point propagator at `sp`.
- `make_volume_src(prop, seed, mix_color, mix_spin, tini)` — volume source.
- `local_sequential_source(res, src, tseq, gammai)` — local sequential source.
- `meson_corr(p0, p1, filename, g0, g1, tini, invmode, info, shift_end, mom)` — meson correlator.
- `corr_dat_create(filename, key_T, dimN, info)` — create a `.dat` correlator file.
- `corr_dat_info(filename, info)` — append info to a `.dat` correlator file.
- `prop4d_conj(prop, rotate)` — complex conjugate a propagator.
- `prop4d_src_gamma(prop, g0, Conj)` — multiply a source gamma matrix.
- `prop4d_sink_gamma(prop, g0, Conj)` — multiply a sink gamma matrix.

## Overview

This module exposes the legacy vector-utility routines that used to live
in the hand-written `qlat/cqlat/vector_utils.cpp` extension.  The
functions operate on the standard `Prop` and `GaugeField` types and
delegate directly to the C++ templates in
`qlat/qlat/include/qlat/vector_utils/`.

It is reached through the ``qlat.c`` re-export module, i.e. the
functions are normally called as ``q.c.<name>(...)`` (aliased as
``qc`` in the examples).

## API Reference

### `diff_gauge(g0, g1)`

Return the maximum difference between the gauge fields `g0` and `g1`.

| Parameter | Type | Description |
|---|---|---|
| `g0` | `GaugeField` | First gauge field. |
| `g1` | `GaugeField` | Second gauge field. |

**Returns:** `float` — the difference measure.

---

### `diff_prop(p0, p1)`

Compare two propagators `p0` and `p1` (no return value; the C++
routine reports the difference).

| Parameter | Type | Description |
|---|---|---|
| `p0` | `Prop` | First propagator. |
| `p1` | `Prop` | Second propagator. |

---

### `load_gwu_link(g0, path)`

Read `g0` from `path` in the gwu link format.

---

### `save_gwu_prop(prop, path)` / `load_gwu_prop(prop, path)`

Write / read a propagator in the gwu format.

---

### `save_gwu_noiP(prop, path)` / `load_gwu_noiP(prop, path)`

Write / read the `noiP` (diagonal noise) projection of a propagator.

---

### `load_qlat_link(g0, path)` / `save_qlat_link(g0, path)`

Read / write a gauge field in the qlat format.

---

### `random_point_src(prop, seed=0)`

Fill `prop` with a random point source.

---

### `make_point_prop(prop, sp=None)`

Fill `prop` with the point propagator located at `sp` (a coordinate
list or `Coordinate`; defaults to the origin).

---

### `make_volume_src(prop, seed=0, mix_color=0, mix_spin=0, tini=-1)`

Fill `prop` with a volume source.  `seed == -1` produces an all-ones
wall source; `mix_color` / `mix_spin` select the colour/spin mixing;
`tini` restricts the source to time slices `>= tini`.

---

### `local_sequential_source(res, src, tseq, gammai=-1)`

Build the local sequential source `res` from `src` for the time
sequence `tseq` and gamma index `gammai`.

---

### `meson_corr(p0, p1, filename, g0, g1, tini=0, invmode=1, info="NONE", shift_end=1, mom=None)`

Compute the meson correlator between `p0` and `p1` and write it to
`filename`.

| Parameter | Type | Description |
|---|---|---|
| `p0`, `p1` | `Prop` | Propagators. |
| `filename` | `str` | Output file name. |
| `g0`, `g1` | `int` | Gamma matrix indices. |
| `tini` | `int` | Initial time slice. |
| `invmode` | `int` | Inversion mode. |
| `info` | `str` | Free-form info string. |
| `shift_end` | `int` | End-shift flag. |
| `mom` | coordinate | Momentum (defaults to zero). |

---

### `corr_dat_create(filename, key_T, dimN, info="NONE")`

Create a `.dat` correlator file with the given key/dimension metadata.

---

### `corr_dat_info(filename, info="NONE")`

Append the info string to an existing `.dat` correlator file.

---

### `prop4d_conj(prop, rotate=1)`

Replace `prop` by its complex conjugate (CPS gamma convention).

---

### `prop4d_src_gamma(prop, g0=0, Conj=0)`

Multiply `prop` by the source gamma matrix `g0` (optionally conjugated).

---

### `prop4d_sink_gamma(prop, g0=0, Conj=0)`

Multiply `prop` by the sink gamma matrix `g0` (optionally conjugated).

## Examples

```python
import qlat as q
import qlat.c as qc

q.begin_with_mpi([[1, 1, 1, 4]])

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))

p0 = q.Prop(geo)
p0.set_rand(q.RngState("p0"))
p1 = q.Prop(geo)
p1 @= p0

qc.save_qlat_link(q.GaugeField(geo), "links.qlat")
qc.prop4d_sink_gamma(p1, 0, 0)
qc.diff_prop(p0, p1)

q.end_with_mpi()
```
