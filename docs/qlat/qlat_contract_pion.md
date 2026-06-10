# `qlat.contract_pion` — Pion Two-Point Contraction

Source: `qlat/qlat/contract_pion.py`

> **Note:** Update this document when updating the source file.

## Outline

- `contract_pion_field(prop, tslice)` — pion correlator on one time slice.

## Overview

This module computes the pion two-point correlation function from a
single quark propagator on a specified lattice time slice.  It supports
both dense `Prop` (Propagator4d) and selected `SelProp` inputs, dispatching
to the appropriate C++ backend in each case.

The result is returned as a `LatData` object containing the correlator
values.  The `@timer` decorator provides automatic profiling.

## API Reference

### `contract_pion_field(prop, tslice)`

Compute the pion two-point function from `prop` on time slice `tslice`.

| Parameter | Type | Description |
|---|---|---|
| `prop` | `Prop` or `SelProp` | Quark propagator (dense or selected). |
| `tslice` | `int` | Lattice time slice index. |

**Returns:** `LatData` — pion correlator for the given time slice.

**Raises:** `Exception` if `prop` is not `Prop` or `SelProp`.

**Decorated with:** `@timer`

## Examples

```python
import qlat as q

q.begin_with_mpi([1, 1, 1, 4])

# Assuming prop is a loaded Prop or SelProp object:
# ld = q.contract_pion_field(prop, tslice=0)
# print(ld)

q.end_with_mpi()
```
