# `qlat_scripts.v1.params_meas` — Physical Measurement Parameters

Source: `qlat/qlat_scripts/v1/params_meas.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Lattice Spacing](#lattice-spacing)
3. [Meson Masses](#meson-masses)
4. [Renormalization Constants (Vector/ Axial)](#renormalization-constants)
5. [Residual Mass](#residual-mass)
6. [Quark Masses](#quark-masses)
7. [Mass Renormalization Factors](#mass-renormalization-factors)
8. [Examples](#examples)

---

## Overview

This module defines physical measurement parameters for each ensemble used in continuum-limit extrapolations and physical-unit conversions. All values are set via `set_param` from `rbc_ukqcd_params` and accessed at runtime via `get_param`.

The parameters include:

- Inverse lattice spacing `a_inv_gev` for converting lattice quantities to GeV
- Meson masses (`m_pi`, `m_kk`) in lattice units
- Renormalization constants (`zz_vv`, `zz_aa`) for vector and axial currents
- Residual mass `m_res` for domain-wall fermions
- Light and heavy quark masses (`m_l`, `m_h`)
- Mass renormalization factors (`zz_m_l`, `zz_m_h`) for converting to MSbar scheme

---

## Lattice Spacing

**Tag:** `a_inv_gev`

Inverse lattice spacing in GeV, used to convert lattice-scale quantities to physical units.

| Ensemble | `a_inv_gev` (GeV) |
|----------|-------------------|
| `test-4nt16` | 1.0 |
| `24D`, `32D`, `48D` | 1.023 |
| `48I` | 1.730 |
| `64I`, `64I-pq` | 2.359 |
| `32IfineH` | 3.148 |

---

## Meson Masses

### `m_pi`

Pion mass in lattice units.

| Ensemble | `m_pi` |
|----------|--------|
| `test-4nt16` | 0.2 |
| `24D`, `32D`, `48D` | 0.13975 |
| `48I` | 0.08049 |
| `32IfineH` | 0.11790 |

### `m_kk`

Kaon mass in lattice units.

| Ensemble | `m_kk` |
|----------|--------|
| `test-4nt16` | 0.4 |
| `24D`, `32D`, `48D` | 0.504154 |
| `48I` | 0.28853 |
| `32IfineH` | 0.17720 |

---

## Renormalization Constants

### `zz_vv`

Vector current renormalization factor.

| Ensemble | `zz_vv` |
|----------|---------|
| `test-4nt16` | 0.7 |
| `24D`, `32D`, `48D` | 0.72672 |
| `48I` | 0.71076 |
| `64I`, `64I-pq` | 0.74293 |
| `32IfineH` | 0.77700 |

### `zz_aa`

Axial current renormalization factor.

| Ensemble | `zz_aa` |
|----------|---------|
| `test-4nt16` | 0.7 |
| `24D`, `32D`, `48D` | 0.73457 |
| `48I` | 0.71191 |
| `64I`, `64I-pq` | 0.74341 |
| `32IfineH` | 0.77779 |

---

## Residual Mass

**Tag:** `m_res`

Domain-wall fermion residual mass in lattice units.

| Ensemble | `m_res` |
|----------|---------|
| `test-4nt16` | 0.001 |
| `24D`, `32D`, `48D` | 0.0022824 |
| `48I` | 0.0006102 |
| `64I`, `64I-pq` | 0.0003116 |
| `32IfineH` | 0.0006296 |

---

## Quark Masses

### `m_l`

Light quark mass in lattice units.

| Ensemble | `m_l` |
|----------|-------|
| `test-4nt8`, `test-4nt16` | 0.01 |
| `24D`, `32D`, `48D` | 0.00107 |
| `48I` | 0.00078 |
| `64I` | 0.000678 |
| `64I-pq` | 0.0006203 |
| `32IfineH` | 0.0047 |

### `m_h`

Heavy (strange) quark mass in lattice units.

| Ensemble | `m_h` |
|----------|-------|
| `test-4nt8`, `test-4nt16` | 0.04 |
| `24D`, `32D`, `48D` | 0.0850 |
| `48I` | 0.0362 |
| `64I` | 0.02661 |
| `64I-pq` | 0.02539 |
| `32IfineH` | 0.0186 |

---

## Mass Renormalization Factors

These convert bare lattice quark masses to the MSbar scheme at 3 GeV. Reference: *Physical Review D* **93**, 074505 (2016).

### `zz_m_l`

Light quark mass renormalization factor.

| Ensemble | `zz_m_l` |
|----------|----------|
| `64I`, `64I-pq` | 2.997 / 2.198 |
| `48I` | 2.997 / 2.198 * 0.9715 |

### `zz_m_h`

Heavy quark mass renormalization factor.

| Ensemble | `zz_m_h` |
|----------|----------|
| `64I`, `64I-pq` | 81.64 / 60.62 |
| `48I` | 81.64 / 60.62 * 0.9628 |

### `zz_ss_l` / `zz_ss_h`

Inverse mass renormalization factors (`1 / zz_m_l` and `1 / zz_m_h`), computed from the above values.

---

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])
import qlat_scripts.v1 as qs

job_tag = "48I"

# Get inverse lattice spacing
a_inv = qs.get_param(job_tag, "a_inv_gev")
print(f"a^{{-1}} = {a_inv} GeV")

# Get pion mass in lattice units
m_pi = qs.get_param(job_tag, "m_pi")
# Convert to MeV
m_pi_mev = m_pi * a_inv * 1000
print(f"m_pi = {m_pi_mev:.1f} MeV")

# Get renormalization constants
zz_vv = qs.get_param(job_tag, "zz_vv")
zz_aa = qs.get_param(job_tag, "zz_aa")
print(f"Z_V = {zz_vv}, Z_A = {zz_aa}")

q.end_with_mpi()
```
