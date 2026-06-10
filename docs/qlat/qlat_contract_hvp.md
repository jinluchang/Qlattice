# `qlat.contract_hvp` — Conserved-Point HVP Contraction (Single Time Slice)

Source: `qlat/qlat/contract_hvp.py`

> **Note:** Update this document when updating the source file.

## Outline

- `contract_chvp3_field(prop1, prop2, tslice)` — contract two propagators on one time slice.

## Overview

This module provides a single high-level wrapper around the C++ routine
`contract_chvp3_sfield`.  It computes the conserved-point hadronic vacuum
polarization (HVP) correlator between two propagators on a specified
lattice time slice and returns the result as a `LatData` object.

The `@timer` decorator records wall-clock time for profiling.

## API Reference

### `contract_chvp3_field(prop1, prop2, tslice)`

Compute the conserved-point HVP on time slice `tslice`.

| Parameter | Type | Description |
|---|---|---|
| `prop1` | `SelProp` / field-like | First propagator (selected or dense). |
| `prop2` | `SelProp` / field-like | Second propagator (selected or dense). |
| `tslice` | `int` | Lattice time slice index. |

**Returns:** `LatData` — the HVP correlator data for the given time slice.

**Decorated with:** `@timer`

## Examples

```python
import qlat as q

q.begin_with_mpi([1, 1, 1, 4])

# Assume prop1 and prop2 are loaded SelProp objects on a given time slice
# ld = q.contract_chvp3_field(prop1, prop2, tslice=3)
# print(ld)

q.end_with_mpi()
```
