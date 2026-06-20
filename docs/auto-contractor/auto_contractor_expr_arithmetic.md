# `auto_contractor.expr_arithmetic` — Symbolic Expression Arithmetic

Source: `qlat/auto_contractor/expr_arithmetic.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Core Classes](#core-classes)
   - [Factor](#factor)
   - [Term](#term)
   - [Expr](#expr)
3. [Interface Functions](#interface-functions)
4. [Simplification](#simplification)
5. [Examples](#examples)

---

## Overview

`expr_arithmetic` provides a symbolic arithmetic layer for manipulating
expression coefficients in the auto-contractor framework. It represents
expressions as sums of `Term` objects, each containing a sympy coefficient
and a list of `Factor` objects (code segments or variables). The module
supports algebraic simplification, Python code compilation, and
variable tracking.

## Core Classes

### Factor

A leaf-level code segment or variable.

| Attribute | Type | Description |
|---|---|---|
| `code` | `str` | The string representation of the factor |
| `variables` | `list[str]` | Variables found in the code |
| `otype` | `str` | `"Var"` if identifier, `"Expr"` otherwise |

```python
f = Factor("a + b")       # otype = "Expr"
f = Factor("x")           # otype = "Var"
```

### Term

A product of `Factor` objects multiplied by a coefficient.

| Attribute | Type | Description |
|---|---|---|
| `coef` | number or sympy expr | The scalar coefficient |
| `factors` | `list[Factor]` | The symbolic factors |

### Expr

A sum of `Term` objects. Supports arithmetic operators (`+`, `-`, `*`)
and simplification.

| Attribute | Type | Description |
|---|---|---|
| `terms` | `list[Term]` | The summands |

## Interface Functions

| Function | Description |
|---|---|
| `mk_expr(x)` | Convert `x` (int, float, str, Factor, Term, Expr) to `Expr` |
| `mk_fac(x)` | Alias for `mk_expr`; creates a factor usable in auto-contractor |
| `mk_sym(x)` | Return `sympy.simplify(x)` |
| `compile_py(x, var_dict=None)` | Compile expression to Python code string |
| `is_zero(x)` | Check if `x` is zero |

## Simplification

| Function | Description |
|---|---|
| `simplified_ea(x)` | Structure-only simplification (sort, combine, drop zeros) |
| `simplified_coef_ea(x)` | Sympy coefficient simplification + structure simplification |

Both return `Expr` or `int` (if the result is a pure scalar).

## Examples

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])

from qlat.auto_contractor.expr_arithmetic import mk_expr, mk_fac, simplified_ea

a = mk_expr(1)
b = mk_expr(2)
c = mk_expr(mk_fac("a + b"))
d = mk_expr(mk_fac("a + b"))
result = c * b + c + d * d
print(result)

result = simplified_ea(result)
print(result)
print(result.compile_py())

q.end_with_mpi()
```
