# `auto_contractor.operators` — Hadronic Operator Construction

Source: `qlat/auto_contractor/operators.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Quark Bilinears](#quark-bilinears)
3. [Meson Operators](#meson-operators)
   - [Pseudoscalar Mesons](#pseudoscalar-mesons)
   - [Scalar Mesons](#scalar-mesons)
   - [Vector Mesons](#vector-mesons)
4. [Two-Meson Operators](#two-meson-operators)
   - [Pi-Pi](#pi-pi)
   - [K-K](#k-k)
   - [K-Pi](#k-pi)
5. [Currents](#currents)
6. [Four-Quark Operators](#four-quark-operators)
7. [Baryon Operators](#baryon-operators)
8. [Examples](#examples)

---

## Overview

`operators` builds symbolic hadronic operators for lattice QCD correlation
functions. It provides factory functions for mesons, baryons, currents, and
four-quark operators (including $\Delta S = 1$ four-quark operators $Q_1$–$Q_{10}$
in both standard and $(8,1)$ representations). All operators are constructed
using `Qb`, `Qv`, `G`, and `Bfield` from the `wick` module and support the
`is_dagger` flag for Hermitian conjugation.

Pi and K operators follow Eq.(103,122) of [Christ et al., Phys. Rev. D 101 (2020) 014506](https://arxiv.org/abs/1908.08640).

## Quark Bilinears

Low-level bilinear constructors used by higher-level operators.

| Function | Definition | Description |
|---|---|---|
| `mk_scalar(f1, f2, p)` | $\bar q_1 q_2(p)$ | Scalar bilinear |
| `mk_scalar5(f1, f2, p)` | $\bar q_1 \gamma_5 q_2(p)$ | Pseudoscalar bilinear |
| `mk_vec_mu(f1, f2, p, mu)` | $\bar q_1 \gamma_\mu q_2(p)$ | Vector bilinear |
| `mk_vec5_mu(f1, f2, p, mu)` | $\bar q_1 \gamma_\mu\gamma_5 q_2(p)$ | Axial-vector bilinear |
| `mk_meson(f1, f2, p)` | $i\bar q_1 \gamma_5 q_2(p)$ | Meson bilinear (with $i$) |

## Meson Operators

### Pseudoscalar Mesons

| Function | Particle | Definition |
|---|---|---|
| `mk_pi_0(p)` | $\pi^0$ | $\frac{i}{\sqrt{2}}(\bar u\gamma_5 u - \bar d\gamma_5 d)$ |
| `mk_pi_p(p)` | $\pi^+$ | $i\bar u\gamma_5 d$ |
| `mk_pi_m(p)` | $\pi^-$ | $-i\bar d\gamma_5 u$ |
| `mk_k_p(p)` | $K^+$ | $i\bar u\gamma_5 s$ |
| `mk_k_m(p)` | $K^-$ | $-i\bar s\gamma_5 u$ |
| `mk_k_0(p)` | $K^0$ | $i\bar d\gamma_5 s$ |
| `mk_k_0_bar(p)` | $\bar K^0$ | $-i\bar s\gamma_5 d$ |
| `mk_eta_l(p)` | $\eta_l$ | $\frac{i}{\sqrt{2}}(\bar u\gamma_5 u + \bar d\gamma_5 d)$ |
| `mk_eta_s(p)` | $\eta_s$ | $i\bar s\gamma_5 s$ |

### Scalar Mesons

| Function | Particle | Definition |
|---|---|---|
| `mk_a0_0(p)` | $a_0^0$ | $\frac{1}{\sqrt{2}}(\bar u u - \bar d d)$ |
| `mk_a0_p(p)` | $a_0^+$ | $\bar u d$ |
| `mk_a0_m(p)` | $a_0^-$ | $\bar d u$ |
| `mk_sigma(p)` | $\sigma$ | $\frac{1}{\sqrt{2}}(\bar u u + \bar d d)$ |
| `mk_kappa_p(p)` | $\kappa^+$ | $\bar u s$ |
| `mk_kappa_m(p)` | $\kappa^-$ | $\bar s u$ |
| `mk_kappa_0(p)` | $\kappa^0$ | $\bar d s$ |
| `mk_kappa_0_bar(p)` | $\bar\kappa^0$ | $\bar s u$ |

### Vector Mesons

| Function | Definition |
|---|---|
| `mk_k_p_star_mu(p, mu)` | $\bar u \gamma_\mu s$ |
| `mk_k_m_star_mu(p, mu)` | $\bar s \gamma_\mu u$ |
| `mk_k_0_star_mu(p, mu)` | $\bar d \gamma_\mu s$ |
| `mk_k_0_star_bar_mu(p, mu)` | $\bar s \gamma_\mu d$ |

## Two-Meson Operators

### Pi-Pi

Isospin decompositions for two-pion operators:

| Function | Isospin | Description |
|---|---|---|
| `mk_pipi_i22(p1, p2)` | $I=2, I_z=+2$ | $\pi^+\pi^+$ |
| `mk_pipi_i21(p1, p2)` | $I=2, I_z=+1$ | $\frac{1}{\sqrt{2}}(\pi^+\pi^0 + \pi^0\pi^+)$ |
| `mk_pipi_i20(p1, p2)` | $I=2, I_z=0$ | $\frac{1}{\sqrt{6}}(2\pi^0\pi^0 + \pi^-\pi^+ + \pi^+\pi^-)$ |
| `mk_pipi_i11(p1, p2)` | $I=1, I_z=+1$ | $\frac{1}{\sqrt{2}}(\pi^+\pi^0 - \pi^0\pi^+)$ |
| `mk_pipi_i10(p1, p2)` | $I=1, I_z=0$ | $\frac{1}{\sqrt{2}}(\pi^+\pi^- - \pi^-\pi^+)$ |
| `mk_pipi_i0(p1, p2)` | $I=0$ | $\frac{1}{\sqrt{3}}(-\pi^0\pi^0 + \pi^-\pi^+ + \pi^+\pi^-)$ |

### K-K

| Function | Isospin | Description |
|---|---|---|
| `mk_kk_i11(p1, p2)` | $I=1, I_z=+1$ | $K^+\bar K^0$ |
| `mk_kk_i10(p1, p2)` | $I=1, I_z=0$ | $\frac{1}{\sqrt{2}}(-K^0\bar K^0 + K^+K^-)$ |
| `mk_kk_i0(p1, p2)` | $I=0$ | $\frac{1}{\sqrt{2}}(K^0\bar K^0 + K^+K^-)$ |
| `mk_k0k0bar(p1, p2)` | — | $K^0\bar K^0$ |

### K-Pi

| Function | Isospin | Description |
|---|---|---|
| `mk_kpi_0_i1half(p1, p2)` | $I=1/2$ | $\frac{1}{\sqrt{3}}K^0\pi^0 + \frac{\sqrt{2}}{\sqrt{3}}K^+\pi^-$ |
| `mk_kpi_p_i1half(p1, p2)` | $I=1/2$ | $\frac{\sqrt{2}}{\sqrt{3}}K^0\pi^+ + \frac{1}{\sqrt{3}}K^+\pi^0$ |
| `mk_kpi_m_i3halves(p1, p2)` | $I=3/2$ | $K^0\pi^-$ |
| `mk_kpi_0_i3halves(p1, p2)` | $I=3/2$ | $-\frac{\sqrt{2}}{\sqrt{3}}K^0\pi^0 + \frac{1}{\sqrt{3}}K^+\pi^-$ |
| `mk_kpi_p1_i3halves(p1, p2)` | $I=3/2$ | $-\frac{1}{\sqrt{3}}K^0\pi^+ + \frac{\sqrt{2}}{\sqrt{3}}K^+\pi^0$ |
| `mk_kpi_p2_i3halves(p1, p2)` | $I=3/2$ | $K^+\pi^+$ |

Two-meson functions accept `is_sym=True` to symmetrize over momentum arguments.

## Currents

| Function | Definition |
|---|---|
| `mk_j5pi_mu(p, mu)` | $\bar d\gamma_\mu\gamma_5 u$ |
| `mk_j5k_mu(p, mu)` | $\bar s\gamma_\mu\gamma_5 u$ |
| `mk_j5km_mu(p, mu)` | $-\bar u\gamma_\mu\gamma_5 s$ |
| `mk_j5eta_l_mu(p, mu)` | $\frac{1}{\sqrt{2}}(\bar u\gamma_\mu\gamma_5 u + \bar d\gamma_\mu\gamma_5 d)$ |
| `mk_j5eta_s_mu(p, mu)` | $\bar s\gamma_\mu\gamma_5 s$ |
| `mk_jpi_mu(p, mu)` | $\bar d\gamma_\mu u$ |
| `mk_jk_mu(p, mu)` | $\bar s\gamma_\mu u$ |
| `mk_j_mu(p, mu)` | $\frac{2}{3}\bar u\gamma_\mu u - \frac{1}{3}\bar d\gamma_\mu d - \frac{1}{3}\bar s\gamma_\mu s$ |
| `mk_jl_mu(p, mu)` | $\frac{2}{3}\bar u\gamma_\mu u - \frac{1}{3}\bar d\gamma_\mu d$ |
| `mk_js_mu(p, mu)` | $-\frac{1}{3}\bar s\gamma_\mu s$ |
| `mk_j0_mu(p, mu)` | $\frac{1}{\sqrt{2}}(\bar u\gamma_\mu u + \bar d\gamma_\mu d)$ (I=0) |
| `mk_j10_mu(p, mu)` | $\frac{1}{\sqrt{2}}(\bar u\gamma_\mu u - \bar d\gamma_\mu d)$ (I=1) |
| `mk_j11_mu(p, mu)` | $\bar u\gamma_\mu d$ (I=1) |
| `mk_j1n1_mu(p, mu)` | $-\bar d\gamma_\mu u$ (I=1) |

## Four-Quark Operators

### Basic Structures

| Function | Structure |
|---|---|
| `mk_4qOp_VV(f1,f2,f3,f4,p)` | $(\bar f_1\gamma_\mu f_2)(\bar f_3\gamma_\mu f_4)$ |
| `mk_4qOp_VA(...)` | $(\bar f_1\gamma_\mu f_2)(\bar f_3\gamma_\mu\gamma_5 f_4)$ |
| `mk_4qOp_AV(...)` | $(\bar f_1\gamma_\mu\gamma_5 f_2)(\bar f_3\gamma_\mu f_4)$ |
| `mk_4qOp_AA(...)` | $(\bar f_1\gamma_\mu\gamma_5 f_2)(\bar f_3\gamma_\mu\gamma_5 f_4)$ |
| `mk_4qOp_SS(...)` | $(\bar f_1 f_2)(\bar f_3 f_4)$ |
| `mk_4qOp_SP(...)` | $(\bar f_1 f_2)(\bar f_3\gamma_5 f_4)$ |
| `mk_4qOp_PS(...)` | $(\bar f_1\gamma_5 f_2)(\bar f_3 f_4)$ |
| `mk_4qOp_PP(...)` | $(\bar f_1\gamma_5 f_2)(\bar f_3\gamma_5 f_4)$ |
| `mk_4qOp_LL(...)` | VV - VA - AV + AA |
| `mk_4qOp_LR(...)` | VV + VA - AV - AA |

### $\Delta S = 1$ Operators ($Q_1$–$Q_{10}$)

| Function | Description |
|---|---|
| `mk_Qsub(p)` | Subtraction operator |
| `mk_Q1(p)` – `mk_Q10(p)` | Standard $\Delta S = 1$ four-quark operators |
| `mk_Q1_b81(p)` – `mk_Q8_b81(p)` | $(8,1)$ representation operators |

## Baryon Operators

| Function | Description |
|---|---|
| `mk_baryon(f1,f2,f3,p,spin,baryon_type)` | Generic baryon from three quarks |
| `mk_proton(p, spin)` | Proton ($uud$) |
| `mk_neutron(p, spin)` | Neutron ($ddu$) |
| `mk_baryon3(f1,f2,f3,p,spin,baryon_type)` | Spin-3/2 baryon |
| `mk_omega(p, spin)` | Omega baryon ($sss$) |

`spin` values: `"u"`, `"d"` for spin-1/2; `"u3"`, `"u1"`, `"d1"`, `"d3"` for spin-3/2.
`baryon_type`: `"std"`, `"pos"` for spin-1/2; `"std3"`, `"pos3"` for spin-3/2.

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])

from qlat.auto_contractor.operators import mk_pi_p, mk_pi_m, mk_k_0, mk_Q1
from qlat.auto_contractor.wick import contract_expr, simplified

# Build a kaon-to-pion weak matrix element operator product
expr = mk_pi_p("x1", is_dagger=True) * mk_Q1("x") * mk_k_0("x2")
c_expr = contract_expr(expr)
c_expr = simplified(c_expr)
print(c_expr.show())

q.end_with_mpi()
```
