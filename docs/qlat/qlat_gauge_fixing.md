# `qlat.gauge_fixing` — Landau and Coulomb Gauge Fixing

Source: `qlat/qlat/gauge_fixing.pyx` (C++ core in
`qlat/qlat/include/qlat/gauge-fixing.h` and `qlat/qlat/lib/gauge-fixing.cpp`)

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Algorithm](#algorithm)
3. [Convergence](#convergence)
4. [Python API](#python-api)
5. [MPI and Performance](#mpi-and-performance)
6. [Example and Cross-Check](#example-and-cross-check)

---

## Overview

`qlat.gauge_fixing` implements **Landau gauge fixing** and **Coulomb gauge
fixing** on a distributed `GaugeField` (`Field<ColorMatrix>`) with a plain
local iterative scheme: every site is updated by the local gauge
transformation that increases the gauge functional.  There is **no Fourier
acceleration and no FFT**.  The site loops live in C++; the Python layer only
calls the Cython wrapper.

## Algorithm

The normalized gauge functional is

$$
F = \frac{1}{V \cdot n_{\mathrm{dir}} \cdot 3} \sum_{x,\mu \in \mathrm{dirs}} \operatorname{Re} \operatorname{Tr} U_\mu(x)
\tag{1}
$$

with all four directions for Landau gauge and the three spatial directions
for Coulomb gauge.  The Coulomb functional splits into one independent
three-dimensional Landau functional per time slice, because a Coulomb gauge
transformation acts inside a single time slice on the spatial links; the
resulting transformation is still applied to all four link directions.

At every site the local staple sum and the local gradient matrix are

$$
K(x) = \sum_{\mu \in \mathrm{dirs}} \left( U_\mu(x) + U_\mu(x-\mu)^\dagger \right)
\tag{2}
$$

$$
D(x) = \frac{K(x)^\dagger - K(x)}{2} - \frac{1}{3} \operatorname{Tr} \frac{K(x)^\dagger - K(x)}{2}
\tag{3}
$$

`D(x)` is anti-hermitian and traceless and is the ascent direction of the
functional.  One Jacobi sweep computes `D(x)` on the whole lattice and then
applies the local gauge transformation to every site at once:

$$
g(x) = \exp(\alpha D(x))
\tag{4}
$$

$$
U_\mu(x) \to g(x) U_\mu(x) g(x+\mu)^\dagger
\tag{5}
$$

The matrix exponential is the truncated Taylor (Horner) recurrence of
`qlat_utils.make_matrix_exp` with order 20.  The step size `alpha` is the
over-relaxation parameter of the plain scheme; it defaults to `0.1`, which is
stable and monotone for the tested lattices.  Larger steps converge faster
but become unstable.

## Convergence

The convergence measure is the RMS lattice divergence

$$
\theta = \sqrt{\frac{1}{V \cdot 3} \sum_x \operatorname{Tr} \left( D(x) D(x)^\dagger \right)}
\tag{6}
$$

which vanishes exactly when the gauge condition `D(x) = 0` holds.  The
iteration stops when `theta <= eps` (default `1e-10`) or after `max_iter`
sweeps (default `10000`).  The returned dictionary contains the full
per-sweep functional and divergence history.

## Python API

All functions take a `qlat.GaugeField` (distributed field) and update it in
place, except the getters:

- `q.fix_gauge(gf, dirs, max_iter=10000, eps=1e-10, alpha=0.1, verbose=False)`
  — fix with an explicit list of link directions.
- `q.fix_landau_gauge(gf, max_iter=10000, eps=1e-10, alpha=0.1, verbose=False)`
  — Landau gauge fixing (all four directions).
- `q.fix_coulomb_gauge(gf, max_iter=10000, eps=1e-10, alpha=0.1, verbose=False)`
  — Coulomb gauge fixing (spatial directions, per time slice).
- `q.get_gauge_functional(gf, dirs)` — normalized gauge functional.
- `q.get_gauge_divergence(gf, dirs)` — RMS divergence `theta`.
- `q.get_unitarity_error(gf)` — maximal `|U U^dagger - 1|` matrix element.

The `fix_*` functions return a dictionary with the keys `gauge`, `converged`,
`n_iter`, `max_iter`, `eps`, `alpha`, `functional`, `divergence`,
`functional_list` and `divergence_list`.

## MPI and Performance

The gauge field is distributed over the `qlat` geometry.  Each sweep only
needs the neighbouring links at the node boundaries: the implementation keeps
an expanded one-hop copy of the gauge field and refreshes it with
`refresh_expanded` (whose `CommPlan` is cached by qlat) once per sweep, and
the same is done for the local transformation.  The functional, the divergence
and the unitarity error are reduced with `field_glb_sum`/`field_glb_max`, so
all nodes stay in sync.  Site loops are in C++ (`qacc_for`); there is no
Python loop over sites.

## Example and Cross-Check

The example `examples-py/gauge-fixing.py` runs on a `4^4` lattice with two
MPI processes, fixes both gauges and checks the convergence, the gauge
condition, unitarity, the monotone increase of the functional and the
invariance of the plaquette.  It is registered in `examples-py/Makefile` and
runs through the examples-py test suite / `nixpkgs/run-one-example-py.py`.

`examples-py/gauge-fixing-cross-check.py` is a manual (non-CI) script that
fixes the same deterministic random field with this implementation and with
the pure-NumPy reference package `qlat-ref`, and compares the converged
functional, the divergence, the unitarity error and the fixed link matrices.
