# `qlat.contract_field` — Site-Level HVP Tensor Contraction

Source: `qlat/qlat/contract_field.py`

> **Note:** Update this document when updating the source file.

## Outline

- `contract_chvp_16(prop1, prop2)` — full 16-component HVP field.

## Overview

This module computes the conserved-point hadronic vacuum polarization
(HVP) tensor with all 16 spin-color polarization combinations
(`mu * 4 + nu`, where `mu, nu ∈ {0,1,2,3}`) at every lattice site.

Unlike `contract_hvp.contract_chvp3_field` which returns a `LatData`
correlator on a single time slice, `contract_chvp_16` returns a full
lattice `Field` with 16 complex components per site.  This is useful
when spatial or point-by-point information is needed.

The underlying C++ function implements:

```
chvp_16(x, mu*4+nu) = tr(g5_herm(prop2(x)) * gamma[mu] * prop1(x) * gamma[nu])
```

where `mu` is the sink polarization and `nu` is the source polarization.

## API Reference

### `contract_chvp_16(prop1, prop2)`

Compute the 16-component HVP tensor field.

| Parameter | Type | Description |
|---|---|---|
| `prop1` | `Propagator4d` | First propagator `prop1_x_y`. |
| `prop2` | `Propagator4d` | Second propagator `prop2_x_y`. |

**Returns:** `Field(ElemTypeComplexD)` — field with 16 complex components
per site, indexed as `mu * 4 + nu`.

**Decorated with:** `@timer`

### Component Index Convention

| Index | mu | nu |
|---|---|---|
| 0 | 0 | 0 |
| 1 | 0 | 1 |
| ... | ... | ... |
| 15 | 3 | 3 |

## Examples

```python
import qlat as q

q.begin_with_mpi([1, 1, 1, 4])

# Assuming prop1 and prop2 are loaded Prop objects:
# chvp = q.contract_chvp_16(prop1, prop2)
# print(chvp.total_elems())  # 16 * number of sites on this node

q.end_with_mpi()
```
