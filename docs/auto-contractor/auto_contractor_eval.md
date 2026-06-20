# `auto_contractor.eval` — Expression Evaluation Engine

Source: `qlat/auto_contractor/eval.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [CCExpr Class](#ccexpr-class)
3. [Compilation and Caching](#compilation-and-caching)
4. [Evaluation](#evaluation)
5. [Benchmarking](#benchmarking)
6. [Random Matrix Generators](#random-matrix-generators)
7. [Utility Functions](#utility-functions)
8. [Meson Build Content](#meson-build-content)
9. [Examples](#examples)

---

## Overview

`eval` is the runtime evaluation module of the auto-contractor framework.  It
takes a compiled `CExpr` (from `compile`), generates Python or Cython code,
compiles it if needed, and provides the `eval_cexpr` entry point for computing
correlator values from propagator lookup functions.  It supports AMA
(all-mode averaging) for multi-precision evaluation.

The typical workflow is:

1. Obtain a `CExpr` via `compile.contract_simplify_compile(*exprs)`.
2. Call `cache_compiled_cexpr(calc_cexpr, path)` to optimise, generate code,
   optionally compile Cython, and return a `CCExpr`.
3. Call `eval_cexpr(ccexpr, positions_dict=..., get_prop=...)` for each
   gauge configuration.

## CCExpr Class

`CCExpr` wraps a compiled and loaded expression module.

```python
ccexpr = CCExpr(cexpr_all, module, base_positions_dict=None, options=None)
```

| Attribute | Type | Description |
|---|---|---|
| `cexpr_all` | dict | Contains `cexpr_original`, `cexpr_optimized`, `code_py` |
| `module` | module | Loaded Python/Cython module with `cexpr_function` |
| `base_positions_dict` | dict | Default position mappings merged into every call |
| `options` | dict | `is_cython`, `is_distillation` flags |
| `cexpr_function_bare` | callable | The raw `cexpr_function` from the module |
| `total_sloppy_flops` | int | FLOP count per sloppy evaluation |
| `expr_names` | list[str] | Names of all expressions |
| `diagram_types` | list | Diagram type definitions |
| `positions` | list[str] | Position variable names |

### Methods

| Method | Description |
|---|---|
| `get_expr_names()` | Return expression name list |
| `cexpr_function(positions_dict, get_prop, is_ama_and_sloppy)` | Evaluate, merging `base_positions_dict` |

## Compilation and Caching

### `cache_compiled_cexpr(calc_cexpr, path, *, is_cython=True, is_distillation=False, base_positions_dict=None)`

The main compilation entry point.  Performs the following steps:

1. Call `calc_cexpr()` to obtain the `CExpr` (cached via pickle).
2. Call `cexpr.optimize()` to run CSE (cached via pickle).
3. Generate Python/Cython code via `cexpr_code_gen_py` (cached via pickle).
4. If Cython: write `meson.build`, run `meson setup` and `meson compile`.
5. Load the compiled module and return a `CCExpr`.

All intermediate results are pickled to `path/` for reuse.  The path is
suffixed with `_cy` or `_py` depending on `is_cython`.  Only node 0 performs
the build; other nodes wait for the pickle file to appear.

### `meson_build_content`

A string containing the Meson build definition for compiling Cython extension
modules.  It links against `qlat-utils`, NumPy, Eigen, OpenMP, and zlib.

## Evaluation

### `eval_cexpr(ccexpr, *, positions_dict, get_prop, is_ama_and_sloppy=False)`

Evaluate a compiled expression.

| Parameter | Type | Description |
|---|---|---|
| `ccexpr` | `CCExpr` | Compiled expression object |
| `positions_dict` | dict | Maps position names to `("tag", coordinate)` tuples |
| `get_prop` | callable | `get_prop(ptype, *args)` returning a Wilson matrix or AMA value |
| `is_ama_and_sloppy` | bool | If True, return `(val_ama, val_sloppy)` tuple |

**Return value:** A 1-D `np.array` of complex values (one per expression), or
a tuple of two such arrays when `is_ama_and_sloppy=True`.

The `get_prop` callback is called as:
- `get_prop("l", pos_snk, pos_src)` — for quark propagators
- `get_prop("U", tag, pos, mu)` — for gauge links

### `get_expr_names(ccexpr)`

Return the list of expression names.  Accepts both `CExpr` and `CCExpr`.

### `get_diagram_type_dict(cexpr)`

Return `diagram_type_dict[diagram_type] = name`.  Accepts both `CExpr` and
`CCExpr`.

## Benchmarking

### `benchmark_eval_cexpr(cexpr, *, benchmark_size=10, benchmark_num=10, benchmark_num_ama=2, benchmark_rng_state=None, base_positions_dict=None)`

Run a correctness and performance benchmark of `eval_cexpr`.  Generates random
positions and propagators, evaluates the expression multiple times, and verifies
that results are consistent across runs and between AMA and sloppy modes.

Returns `(check, check_ama)` — lists of complex fingerprints derived from
dotting the result arrays with random check vectors.

## Random Matrix Generators

These functions create random test matrices for benchmarking:

| Function | Return Type (normal) | Return Type (distillation) |
|---|---|---|
| `make_rand_spin_color_matrix(rng_state)` | `WilsonMatrix` (12×12 complex) | `np.ndarray` (4, nc, 4, nc) |
| `make_rand_spin_matrix(rng_state)` | `SpinMatrix` (4×4 complex) | `np.ndarray` (4, 4) |
| `make_rand_color_matrix(rng_state)` | `ColorMatrix` (3×3 complex) | `np.ndarray` (nc, nc) |

When `is_distillation=True`, colour dimension `nc=10` and spin dimension `ns=4`.

## Utility Functions

| Function | Description |
|---|---|
| `benchmark_show_check(check)` | Format a check vector as a space-separated scientific notation string |
| `sqr_component(x)` | Square real and imaginary parts independently |
| `sqrt_component(x)` | Square root of real and imaginary parts independently |
| `sqr_component_array(arr)` | Element-wise `sqr_component` on an array |
| `sqrt_component_array(arr)` | Element-wise `sqrt_component` on an array |

## Examples

```python
import qlat as q
import auto_contractor as ac

q.begin_with_mpi([[1, 1, 1, 4]])

# Build and compile a simple expression
def calc_cexpr():
    expr = ac.Qb("u", "x", "s1", "c1") * ac.G(5, "s1", "s2") * ac.Qv("u", "x", "s2", "c1")
    return ac.contract_simplify_compile(expr, is_isospin_symmetric_limit=True)

ccexpr = ac.cache_compiled_cexpr(calc_cexpr, "cache/pi_test", is_cython=False)

# Evaluate with a dummy propagator
size = q.Coordinate([4, 4, 4, 8])
positions_dict = {
    "size": size,
    "x": ("point", q.Coordinate([0, 0, 0, 0])),
}

def get_prop(ptype, *args):
    if ptype == "U":
        return q.ColorMatrix()
    wm = q.WilsonMatrix()
    arr = wm[:]
    arr[:] = 0
    return wm

result = ac.eval_cexpr(ccexpr, positions_dict=positions_dict, get_prop=get_prop)
print(f"Result: {result}")

q.end_with_mpi()
```
