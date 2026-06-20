# `auto_contractor.wick` — Wick Contraction Engine

Source: `qlat/auto_contractor/wick.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Operator Classes](#operator-classes)
   - [Quark Fields](#quark-fields)
   - [Propagators and Matrices](#propagators-and-matrices)
   - [Traces and Chains](#traces-and-chains)
   - [Baryon Fields](#baryon-fields)
3. [Expression System](#expression-system)
4. [Contraction](#contraction)
5. [Simplification](#simplification)
6. [Utility Functions](#utility-functions)
7. [Examples](#examples)

---

## Overview

`wick` implements the symbolic Wick contraction engine for the auto-contractor
framework. Given products of quark creation/annihilation operators, gamma
matrices, and gauge links, it performs Wick contractions to produce expressions
in terms of quark propagator traces. The result is a sum of terms, each
containing traces (`Tr`) or open chains (`Chain`) of propagators, gamma
matrices, and gauge links, multiplied by symbolic coefficients.

The module also handles baryon operators via the `Bfield` and `BS` (baryon
propagator) classes, supporting spin-1/2 and spin-3/2 baryons with both
standard and positive-parity interpolating fields.

## Operator Classes

### Quark Fields

All inherit from `Qfield` (which inherits from `Op`).

| Class | Otype | Description |
|---|---|---|
| `Qv(f, p, s, c)` | `"Qv"` | Quark annihilation operator $\psi(p)$; acts as $\partial/\partial\bar H$ |
| `Qb(f, p, s, c)` | `"Qb"` | Quark creation operator $\bar\psi(p)$; acts as $\partial/\partial H$ |
| `Hv(f, p, s, c)` | `"Hv"` | Hadronic annihilation operator |
| `Hb(f, p, s, c)` | `"Hb"` | Hadronic creation operator |
| `SHv(f, p, s, c)` | `"SHv"` | Smeared hadronic annihilation |
| `HbS(f, p, s, c)` | `"HbS"` | Smeared hadronic creation |

Parameters: `f` = flavor, `p` = position label, `s` = spin index, `c` = color index.

### Propagators and Matrices

| Class | Otype | Description |
|---|---|---|
| `S(f, p1, p2, ...)` | `"S"` | Quark propagator from `p2` to `p1` with flavor `f` |
| `G(tag, s1, s2)` | `"G"` | Gamma matrix; `tag` = 0,1,2,3 (spatial/temporal) or 5 ($\gamma_5$) |
| `U(tag, p, mu, c1, c2)` | `"U"` | Gauge link at position `p` in direction `mu` |

### Traces and Chains

| Class | Otype | Description |
|---|---|---|
| `Tr(ops, tag)` | `"Tr"` | Trace over spin/color of a product of S, G, U operators |
| `Chain(ops, tag, ...)` | `"Chain"` | Open (non-loop) product of S, G, U operators |

The `tag` attribute indicates what indices are traced: `"sc"` (spin+color), `"s"` (spin only), `"c"` (color only), or `""` (neither).

### Baryon Fields

| Class | Description |
|---|---|
| `Bfield(tag, s1, s2, s3, c1, c2, c3)` | Baryon creation/annihilation tensor |
| `BfieldCoef(st_list_code)` | Spin-tensor coefficients for baryon operators |
| `BS(elem_list, chain_list)` | Baryon propagator (three chained quark propagators) |

Predefined baryon tags in `bfield_tag_dict`:

| Tag | Description |
|---|---|
| `"std-u"`, `"std-d"` | Standard spin-1/2 proton/neutron |
| `"pos-u"`, `"pos-d"` | Positive-parity spin-1/2 |
| `"std3-u3"`, `"std3-u1"`, `"std3-d1"`, `"std3-d3"` | Standard spin-3/2 |
| `"pos3-u3"`, `"pos3-u1"`, `"pos3-d1"`, `"pos3-d3"` | Positive-parity spin-3/2 |

## Expression System

| Class | Description |
|---|---|
| `Term(c_ops, a_ops, coef)` | A single term: coefficient times commutable and anti-commutable operators |
| `Expr(terms, description)` | A sum of terms with optional description label |

Key `Expr` methods:

| Method | Description |
|---|---|
| `sort()` | Sort operators within each term |
| `simplify()` | Full simplification pipeline |
| `collect_traces()` | Identify and form `Tr` and `Chain` objects |
| `combine_terms()` | Merge terms with identical operator structures |
| `simplify_coef()` | Simplify coefficients via sympy |
| `rescale_bs_term()` | Normalize baryon propagator coefficients |
| `round(ndigit)` | Numerically evaluate coefficients to `ndigit` digits |

Arithmetic operators (`+`, `-`, `*`) are supported. Adding a `str` sets the description.

## Contraction

| Function | Description |
|---|---|
| `contract_expr(expr)` | Perform Wick contractions on `Qv`/`Qb` operators |
| `contract_term(term)` | Contract a single term |

The contraction pushes `Qv` and `Qb` operators through the term, replacing
them with propagators `S` via:
- $Qv \cdot HbS \to S$ (with sign tracking for fermion ordering)
- $Qb \cdot SHv \to S$

## Simplification

| Function | Description |
|---|---|
| `simplified(expr)` | Full simplification (isospin limit, traces, combine, scale, simplify coefficients) |
| `collect_traces(ops)` | Identify `Chain` and `Tr` sub-structures from a flat operator list |
| `collect_baryon_props(op_list)` | Identify `BS` (baryon propagator) structures |

## Utility Functions

| Function | Description |
|---|---|
| `mk_fac(x)` | Create a factor expression (delegates to `expr_arithmetic.mk_fac`) |
| `mk_expr(x)` | Convert an `Op`, `Term`, or scalar to `Expr` |
| `S_l(p1, p2)`, `S_s(p1, p2)`, `S_c(p1, p2)` | Propagator shortcuts for light, strange, charm |
| `tr(expr)` | Wrap commutable operators in a `Tr` |
| `gamma(tag)` | Create a gamma matrix expression |
| `gamma_va(tag)` | Vector (`0`–`3`) or axial-vector (`4`–`7`) gamma |
| `gamma_x`, `gamma_y`, `gamma_z`, `gamma_t`, `gamma_5` | Pre-built gamma expressions |

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])

from qlat.auto_contractor.wick import (
    Qb, Qv, G, S, Tr, Expr, Term,
    contract_expr, simplified, mk_expr,
)

# Build a pion correlator: tr[g5 S_d(x,y) g5 S_u(y,x)]
expr = (
    1
    * Qb("d", "x1", "s1", "c1")
    * G(5, "s1", "s2")
    * Qv("u", "x1", "s2", "c1")
    * Qb("u", "x2", "s3", "c2")
    * G(5, "s3", "s4")
    * Qv("d", "x2", "s4", "c2")
)

# Perform Wick contractions
c_expr = contract_expr(expr)
c_expr = simplified(c_expr, is_isospin_symmetric_limit=False)
print(c_expr.show())

q.end_with_mpi()
```
