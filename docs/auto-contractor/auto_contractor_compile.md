# `auto_contractor.compile` — Expression Compiler

Source: `qlat/auto_contractor/compile.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Key Classes](#key-classes)
   - [Var](#var)
   - [CExpr](#cexpr)
   - [CExprCodeGenPy](#cexprcodegenpy)
3. [Interface Functions](#interface-functions)
4. [Common Sub-Expression Elimination](#common-sub-expression-elimination)
5. [Diagram Type Filtering](#diagram-type-filtering)
6. [Code Generation](#code-generation)
7. [Helper Functions](#helper-functions)
8. [Examples](#examples)

---

## Overview

`compile` is the central compilation module of the auto-contractor framework.
It takes Wick-contracted `Expr` objects (produced by `wick.contract_expr`) and
transforms them into an optimised intermediate representation (`CExpr`) through
common sub-expression elimination.  The `CExpr` can then be serialised to
Python or Cython source code that evaluates the contractions efficiently at
runtime.

The typical workflow is:

1. Build symbolic expressions with quark-field operators (`Qb`, `Qv`, `S`, `G`, …).
2. Call `contract_simplify_compile(*exprs)` to Wick-contract, simplify, and
   compile into a `CExpr`.
3. Pass the `CExpr` to `eval.cache_compiled_cexpr` which calls `cexpr.optimize()`,
   generates code, and produces a loadable `CCExpr`.

## Key Classes

### `Var`

A variable placeholder node in the compiled expression tree.

```python
class Var(Op):
    def __init__(self, name: str)
```

`Var` objects replace repeated sub-expressions (propagators, colour matrices,
traces, chains, products) with named variables during compilation.

### `CExpr`

The compiled expression container.  Created by `mk_cexpr` and populated by
`CExpr.optimize()` / `CExpr.collect_op()`.

| Attribute | Type | Description |
|---|---|---|
| `diagram_types` | `list[tuple]` | Named diagram types |
| `positions` | `list[str]` | Position variable names |
| `variables_factor_intermediate` | `list[tuple]` | Intermediate factor definitions |
| `variables_factor` | `list[tuple]` | Final factor (coefficient) definitions |
| `variables_prop` | `list[tuple]` | Propagator variable definitions |
| `variables_color_matrix` | `list[tuple]` | Colour matrix variable definitions |
| `variables_prod` | `list[tuple]` | Common product sub-expression definitions |
| `variables_chain` | `list[tuple]` | Chain variable definitions |
| `variables_tr` | `list[tuple]` | Trace variable definitions |
| `variables_baryon_prop` | `list[tuple]` | Baryon propagator definitions |
| `named_terms` | `list[tuple]` | Named term definitions |
| `named_exprs` | `list[tuple]` | Named expression definitions |

Key methods:

| Method | Description |
|---|---|
| `optimize()` | Run full CSE pipeline (must be called before evaluation) |
| `collect_op()` | Perform common sub-expression elimination |
| `copy()` | Return a deep copy |
| `get_expr_names()` | Return list of expression names |

### `CExprCodeGenPy`

Code generator that emits Python or Cython source from an optimised `CExpr`.

```python
gen = CExprCodeGenPy(cexpr, is_cython=True, is_distillation=False)
code = gen.code_gen()
```

The generated code defines `cexpr_function(positions_dict, get_prop, ...)` and
its helper functions.  FLOP counts are tracked per matrix operation type
(13536 for sc×sc, 4320 for sc×s, 480 for s×s, etc.).

## Interface Functions

### `mk_cexpr(*exprs, diagram_type_dict=None)`

Build a `CExpr` from already-contracted expressions.  Assigns diagram type
names, term names, and collects positions.

### `contract_simplify(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=None)`

Wick-contract and simplify each expression.  Accepts raw `Expr` objects or
`(expr, *included_types)` tuples for diagram type filtering.

### `compile_expr(*exprs, diagram_type_dict=None)`

Wrap contracted/simplified expressions into a `CExpr`.

### `contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=None)`

Convenience function that calls `contract_simplify` then `compile_expr`.  This
is the primary entry point; its result can be passed directly to
`cache_compiled_cexpr`.

### `display_cexpr(cexpr)`

Return a human-readable string representation of a `CExpr` showing all
variable definitions, terms, and expressions.

### `cexpr_code_gen_py(cexpr, *, is_cython=True, is_distillation=False)`

Generate Python or Cython source code for a `CExpr`.  If `is_distillation` is
True, `is_cython` must be False.

## Common Sub-Expression Elimination

The `collect_op` pipeline performs several CSE passes:

1. **Factor collection** (`collect_and_optimize_factor_in_cexpr`) — Extract
   numerical coefficients into named variables and find common products/sums.
2. **Propagator collection** (`collect_prop_in_cexpr`) — Deduplicate `S` (propagator)
   operators into `V_S_*` variables.
3. **Colour matrix collection** (`collect_color_matrix_in_cexpr`) — Deduplicate
   `U` operators into `V_U_*` variables.
4. **Chain collection** (`collect_chain_in_cexpr`) — Deduplicate `Chain`
   operators into `V_chain_*` variables.
5. **Trace collection** (`collect_tr_in_cexpr`) — Deduplicate `Tr` operators
   into `V_tr_*` variables.
6. **Baryon propagator collection** (`collect_baryon_prop_in_cexpr`) — Deduplicate
   `BS` operators into `V_bs_*` variables.
7. **Product sub-expression collection** (`collect_subexpr_in_cexpr`) — Find
   common adjacent operator pairs inside traces and replace them with `V_prod_*`
   variables.

## Diagram Type Filtering

`filter_diagram_type(expr, diagram_type_dict, included_types)` partitions an
expression by its diagram topology.  Each term's diagram type is computed from
its propagator connectivity (with permutation symmetry).  Terms whose type maps
to `None` in `diagram_type_dict` are dropped.

The `included_types` parameter controls which named types appear in each output
expression.  Each entry can be `None` (keep all), a single type name string,
or a list of type name strings.

## Code Generation

`CExprCodeGenPy.code_gen()` emits a self-contained Python (or Cython) module
with:

- `cexpr_function(positions_dict, get_prop, is_ama_and_sloppy)` — top-level entry
- `cexpr_function_get_prop(...)` — fetches propagators and colour matrices
- `cexpr_function_eval(...)` — applies AMA (all-mode averaging) to the sloppy evaluation
- `cexpr_function_eval_with_props(...)` — the core computation kernel
- `cexpr_function_bs_eval_*` — baryon contraction kernels
- `total_sloppy_flops` — estimated FLOP count per sloppy call

## Helper Functions

| Function | Description |
|---|---|
| `get_positions(term)` | Return sorted list of position strings in a term |
| `get_term_diagram_type_info(term)` | Compute the diagram type of a term |
| `show_variable_value(value)` | Human-readable string for a variable's RHS |
| `get_var_name_type(name)` | Return type tag (`V_S`, `V_G`, `V_U`, `V_a`) from a variable name |
| `get_op_type(op)` | Return type string for an `Op` node |

## Examples

```python
import qlat as q
import auto_contractor as ac

q.begin_with_mpi([[1, 1, 1, 4]])

# Build a simple pion correlator expression
expr = ac.Qb("u", "x", "s1", "c1") * ac.G(5, "s1", "s2") * ac.Qv("u", "x", "s2", "c1")

# Contract, simplify, and compile
cexpr = ac.contract_simplify_compile(expr, is_isospin_symmetric_limit=True)
print(ac.display_cexpr(cexpr))

q.end_with_mpi()
```
