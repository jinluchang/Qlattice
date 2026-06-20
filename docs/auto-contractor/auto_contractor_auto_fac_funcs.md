# `auto_contractor.auto_fac_funcs` — Modular Arithmetic Helpers

Source: `qlat/auto_contractor/auto_fac_funcs.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Exported Functions](#exported-functions)
3. [Function List](#function-list)
4. [Examples](#examples)

---

## Overview

`auto_fac_funcs` is a thin re-export module that makes the modular arithmetic
helper functions from `qlat_utils` available under the `auto_contractor`
namespace.  These functions are referenced by name in compiled expression code
(through `auto_fac_funcs_list`) and resolved at evaluation time.

## Exported Functions

All three functions are imported directly from `qlat_utils`:

| Function | Description |
|---|---|
| `rel_mod(x, size)` | Compute `x mod size` (non-negative remainder) |
| `rel_mod_sym(x, size)` | Compute `x mod size` centred around zero (symmetric modulo) |
| `c_rel_mod_sqr(x, size)` | Squared modulus of the relative modular distance |

## Function List

The module exposes `auto_fac_funcs_list`, a plain list of the function names
above.  The compiled expression generator (`compile.cexpr_code_gen_py`) uses
this list to recognise which position variables should be resolved to these
helper functions rather than looked up from `positions_dict`.

```python
auto_fac_funcs_list = [
    "rel_mod",
    "rel_mod_sym",
    "c_rel_mod_sqr",
]
```

## Examples

```python
import qlat as q
import auto_contractor as ac

q.begin_with_mpi([[1, 1, 1, 4]])

from auto_contractor.auto_fac_funcs import rel_mod, rel_mod_sym, c_rel_mod_sqr

print(rel_mod(7, 4))        # 3
print(rel_mod_sym(7, 4))    # -1
print(c_rel_mod_sqr(7, 4))  # 9.0

q.end_with_mpi()
```
