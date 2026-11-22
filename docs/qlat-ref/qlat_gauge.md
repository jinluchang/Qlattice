# `qlat_ref.gauge` — Landau and Coulomb Gauge Fixing (NumPy Reference)

Source: `qlat-ref/qlat_ref/gauge.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Array Layout and Conventions](#array-layout-and-conventions)
3. [Algorithm](#algorithm)
4. [Convergence and Monitoring](#convergence-and-monitoring)
5. [Public Interface](#public-interface)
6. [Examples](#examples)
7. [Self-Checks](#self-checks)

---

## Overview

`qlat_ref.gauge` is a deliberately small, pure-NumPy reference implementation
of **Landau gauge fixing** and **Coulomb gauge fixing** on a
four-dimensional lattice.  It is the reference used to cross-check the
efficient C++/Cython implementation in `qlat`; it is not a performance target.

The iteration scheme is **plain local steepest descent**: every site is
updated at once by the local gauge transformation that increases the gauge
functional, without Fourier acceleration and without FFT.  The whole lattice is
treated with NumPy array operations, so there is no loop over sites in the
code.

## Array Layout and Conventions

A gauge field is a NumPy array `U` with shape `(4, Nt, Nx, Ny, Nz, 3, 3)` and
`dtype` `numpy.complex128`:

    U[mu, t, x, y, z, i, j] = U_mu(t, x, y, z)[i, j]

`U[mu]` is the field of links in direction `mu` (with `0 = t`, `1 = x`,
`2 = y`, `3 = z`), so the site axes of `U[mu]` are the time and the three
spatial directions, and the last two axes are the (row, column) color indices
of the `SU(3)` matrix.

A gauge transformation is a site field of shape `(Nt, Nx, Ny, Nz, 3, 3)`.
The lattice is periodic in every direction: a shift by one hop in direction
`mu` is `numpy.roll(..., 1, axis=mu)`.  For a link field `U[mu]` the
coordinate axis of direction `mu` is `axis=mu`, because `U[mu]` drops the
leading direction axis of `U`.

The gauge transformation convention is

$$
U_\mu(x) \to g(x) U_\mu(x) g(x+\mu)^\dagger
\tag{1}
$$

which is implemented by `apply_gauge_transform`.

## Algorithm

The Landau gauge functional is

$$
F = \frac{1}{V \cdot 4 \cdot 3} \sum_{x,\mu} \operatorname{Re} \operatorname{Tr} U_\mu(x)
\tag{2}
$$

and the Coulomb gauge functional is the same sum restricted to the spatial
directions `i = 1, 2, 3`.  The Coulomb functional splits into one independent
three-dimensional Landau functional per time slice, because a Coulomb gauge
transformation `g(t, x)` acts inside a single time slice on the spatial links.

At every site the local gradient matrix is

$$
D(x) = \frac{K(x)^\dagger - K(x)}{2} - \frac{1}{3} \operatorname{Tr}\frac{K(x)^\dagger - K(x)}{2}
\tag{3}
$$

with the staple sum

$$
K(x) = \sum_{\mu} \left( U_\mu(x) + U_\mu(x-\mu)^\dagger \right)
\tag{4}
$$

the sum running over all four directions for Landau gauge and over the
spatial directions for Coulomb gauge.  `D(x)` is anti-hermitian and traceless,
so the local gauge transformation

$$
g(x) = \exp(\alpha D(x))
\tag{5}
$$

is in `SU(3)` (the matrix exponential is a truncated Taylor series, see
`get_su3_exp`).  The step size `alpha` defaults to `0.1`; it is the
over-relaxation parameter of the plain scheme.  Larger steps converge faster
but become unstable, so `alpha = 0.1` is the safe default.

One sweep computes `D(x)` on the whole lattice and then applies the
transformation to every site at once (Jacobi sweep, no site loop and no
ordering dependence).  The gauge field `U` passed to the fixing functions is
not mutated; a gauge-transformed copy is returned.

## Convergence and Monitoring

The convergence measure is the RMS lattice divergence

$$
\theta = \sqrt{\frac{1}{V \cdot 3} \sum_x \operatorname{Tr} \left( D(x) D(x)^\dagger \right)}
\tag{6}
$$

which vanishes exactly when the gauge condition `D(x) = 0` holds.  The
iteration stops when `theta <= eps` (default `eps = 1e-10`) or after
`max_iter` sweeps (default `10000`).  Both the functional and the divergence
of every sweep are recorded in the returned `info` dictionary, which makes it
easy to monitor the convergence.

## Public Interface

- `get_unit_gauge_field(total_site)` — unit gauge field.
- `get_random_gauge_field(total_site, *, seed=0, scale=1.0)` — deterministic
  random unitary gauge field.
- `get_random_gauge_transform(total_site, *, seed=0, scale=0.3)` — deterministic
  random `SU(3)` gauge transformation.
- `apply_gauge_transform(U, g)` — apply a gauge transformation.
- `get_gauge_functional(U, *, dirs=None)` — normalized gauge functional.
- `get_coulomb_functional_list(U)` — Coulomb functional of each time slice.
- `get_gauge_divergence(U, *, dirs=None)` — RMS divergence `theta`.
- `get_coulomb_divergence_list(U)` — RMS divergence of each time slice.
- `get_local_gauge_gradient(U, *, dirs=None)` — local gradient field `D(x)`.
- `get_staple_sum(U, *, dirs=None)` — staple sum `K(x)`.
- `get_unitarity_error(U)` / `check_unitarity(U, *, eps=1e-12)` — unitarity
  monitoring.
- `fix_landau_gauge(U, *, max_iter=10000, eps=1e-10, alpha=0.1,
  verbose=False)` — Landau gauge fixing.
- `fix_coulomb_gauge(U, *, max_iter=10000, eps=1e-10, alpha=0.1,
  verbose=False)` — Coulomb gauge fixing.

`fix_landau_gauge` and `fix_coulomb_gauge` return a tuple `(U_fixed, info)`.
`info` contains the gauge name, the direction set, the final functional and
divergence, the convergence flag, the number of sweeps and the per-sweep
history.

## Examples

```python
import numpy as np
import qlat_ref as qr

total_site = (2, 2, 2, 2)
U = qr.get_random_gauge_field(total_site, seed=0)
U_fixed, info = qr.fix_landau_gauge(U)
assert info["converged"]
assert qr.get_gauge_divergence(U_fixed) <= info["eps"]
assert qr.check_unitarity(U_fixed, eps=1e-10)

U_coulomb, info_coulomb = qr.fix_coulomb_gauge(U)
assert info_coulomb["converged"]
assert max(qr.get_coulomb_divergence_list(U_coulomb)) <= np.sqrt(U.shape[1]) * info_coulomb["eps"]
```

## Self-Checks

The module owns a runnable self-test:

```bash
cd qlat-ref
python3 -m qlat_ref.tests
```

It checks that the random fields are unitary, that the functional increases
monotonically during the plain iteration, that both gauges converge and
satisfy their gauge condition, that a `1^4` lattice is trivially gauge fixed
(the functional is invariant and the divergence vanishes), that the iteration
is exactly covariant under a global gauge transformation, and that the
converged functional does not depend on the step size.
