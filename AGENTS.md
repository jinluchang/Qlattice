# AGENTS.md — Guide for Coding Agents

## Project Overview

Qlattice is a lattice QCD simulation library. The monorepo contains four Python/C++ packages:
- `qlat-utils` — core utilities (C++17, Cython, Python)
- `qlat` — main lattice library (depends on qlat-utils)
- `qlat-grid` — Grid framework interface
- `qlat-cps` — CPS framework interface

Primary languages: **C++17**, **Cython**, **Python**, **Bash**.

> **Note**: The `applications/` directory contains old code that may be outdated or no longer correct. Do NOT reference it as a source of truth or as an example to follow.
>
> **Note**: The `examples-*` directories (`examples-py/`, `examples-py-gpt/`, `examples-py-cps/`, `examples-cpp/`, `examples-cpp-grid/`) ARE the canonical reference for current API usage and are safe to reference. Keep them up to date when APIs change, and run their tests routinely to catch regressions.

## Build System

Build system: **Nix** which calls **Meson** to build the packages (via `meson-python` for pip packaging). No CMake or npm.

## Build and Testing

### Build qlat with nix

Use `nixpkgs/install-py-local-kernel-with-nix.py` to build qlat via nix. It creates a `./result-py-local` symlink to the nix store path. Build variants: `--variant cuda`, `--variant cudasupport`, `--variant cu`, `--variant clang`, `--variant pypi`, or `--all-variants` for every variant. See `--help` for details.

```bash
./nixpkgs/install-py-local-kernel-with-nix.py
```

### Run a single test

**REQUIREMENT**: All single tests MUST be run via the `nixpkgs/run-one-example*.py` scripts. Do NOT run tests manually or via `make` directly — these scripts set up the nix environment correctly and delegate to `make`:

```bash
./nixpkgs/run-one-example-py.py utils                   # Python test
./nixpkgs/run-one-example-py-gpt.py qtopo-measure-dwf   # GPT test
./nixpkgs/run-one-example-cpp.py simple-1               # C++ test
./nixpkgs/run-one-example-cpp-grid.py grid-with-qlat    # Grid test
./nixpkgs/run-one-example-py-cps.py cps-qlat-io-test    # CPS test
```

See `--help` for build variant options (`--cuda`, `--cudasupport`, `--cu`, `--clang`, `--pypi`).

**Build requirement**: Before running a test, ensure qlat is built with nix. If `./result-py-local` does not exist, or if source code has changed since the last build (check `git status` or file timestamps), build/rebuild first:
```bash
./nixpkgs/install-py-local-kernel-with-nix.py
```

Each script copies sources into `./tmp/examples-*/` and runs the test there. After a run, check `./tmp/examples-*/<name>.py.p/` (Python/GPT/CPS) or `./tmp/examples-*/<name>/build/` (C++/Grid) for log files (`log.full.txt`, `log.txt`, `log.check.txt`).

Tests use **log-comparison**: each test prints `CHECK:` lines compared against reference `.log` files. See each `examples-*/Makefile` for the full test-running workflow.

### Run all tests

**REQUIREMENT**: All tests MUST be run via the `nixpkgs/build-many-qlat-pkgs.py` script below. Do NOT loop through individual tests manually — the `tests` package set builds the `qlat-tests` package for the default version and name, which runs the whole test suite (Python, GPT, CPS, C++ and Grid examples) and fails when any `log.check.txt` differs:

```bash
./nixpkgs/build-many-qlat-pkgs.py --group tests
```

The script sets up the nix and nom caches (falling back to `<repo>/tmp/` when the user directories are not writable) and writes the test package to `$HOME/qlat-build/nix/tests/result` (`<repo>/tmp/qlat-build/nix/tests/result` when `$HOME` is not writable). `--group` may be repeated and combined with `-j`/`--cores`, and `./nixpkgs/build-many-qlat-pkgs.py --list-groups` shows every package set.

**DO NOT** use shell loops like `for test in ... ; do ./nixpkgs/run-one-example-py.py $test ; done` — this is incorrect and bypasses the proper test orchestration.

### Run a new program with the nix-built environment

Always run new programs under `<project_root>/tmp/` (or a sub-directory) to keep test artifacts out of the source tree.

Source `setenv-qlat.sh` from the nix result to configure PATH, PYTHONPATH, LD_LIBRARY_PATH, etc.:

```bash
source result-py-local/bin/setenv-qlat.sh
mkdir -p tmp
cd tmp
mpiexec -n 2 --oversubscribe --bind-to none bash bind-gpu-qlat.sh python3 -m mpi4py ./your-script.py
```

See the `%.log` rule in `examples-py/Makefile` for the full invocation pattern (mpi options, CHECK-line extraction, etc.).

## Code Style — C++

Configuration: `.clang-format` (Google base, Linux braces, left pointers).

### Formatting
- **Indentation**: 2 spaces (no tabs). Vim modeline: `ts=2 sw=2 expandtab`
- **Braces**: Linux style — opening brace on next line for functions, same line for control flow
- **Pointers**: Left-aligned (`int* p`, not `int *p`)
- **Standard**: C++17 (`cpp_std=c++17` in meson.build)

### Naming
- Classes/structs: `PascalCase` — e.g., `GeometryNode`, `Field`, `SelectedField`
- Functions: `snake_case` — e.g., `get_elem`, `coordinate_from_index`
- Variables: `snake_case`; private members end with trailing underscore
- Constants: `UPPER_CASE` or `PascalCase` for enum values
- Template params: `PascalCase` — e.g., `template <class M>`
- Type aliases: `using` with `PascalCase` — e.g., `using GaugeField = GaugeFieldT<>;`

### Includes & Headers
- `#pragma once` for include guards (no `#ifndef` guards)
- Project headers first, then system headers
- All code in `namespace qlat { ... }`

### Error Handling
- Use project macros: `qassert(cond)`, `Qassert(cond)`, `qerr(msg)`, `Qassert_info(cond, {...})`
- Warnings: `warn(msg)` or `displayln_info(ssprintf("WARNING: ..."))`
- Test markers: `displayln_info("CHECK: ...")` for test verification

### Key Macros/Decorators
- `qacc` — accelerator decorator (CPU/GPU portability)
- `TIMER("name")` — scope-based timer
- `API` — DLL export macro

## Code Style — Python

### Imports
```python
import qlat_gpt as qg   # MUST be first if used (initializes GPT environment)
import numpy as np
import qlat as q
import qlat_utils as qu   # if needed
```
`qlat_gpt` must be imported before all other packages when it is used, as it initializes the GPT/Grid runtime environment. After that, standard library, then `qlat` (aliased as `q`). Avoid `import *`; use explicit imports instead.

### Conventions
- Shebang: `#!/usr/bin/env python3` for executable scripts
- `snake_case` for functions/variables, `PascalCase` for classes
- f-strings for formatting
- No type hints in example scripts (but welcomed in library code)

### Test Script Pattern
```python
# Global function definitions BEFORE q.begin_with_mpi()
def my_helper():
    ...

q.begin_with_mpi(size_node_list)
q.json_results_append("test description")
# ... test logic with assert for correctness ...
q.json_results_append("plaq", plaq_value, 1e-10)  # floating-point result
q.json_results_append(f"n_marks = {n_marks}")     # int/bool -> name only
```
- **Functions first**: Define all helper functions before `q.begin_with_mpi()`.
- **Results**: `q.json_results_append(name, value[, eps])` is only for floating-point results (`float`, `complex`, a numeric `numpy.ndarray`). For any non-float result (boolean/flag, integer count, string, exception outcome, plain marker) pass a single string and encode the outcome in it, e.g. `q.json_results_append(f"match = {ok}")`. Never fabricate a float (`float(ok)`, `float(len(x))`, `1.0` for `True`).
- **No intermediate CHECK lines**: Use `q.json_results_append` to record test results. Do NOT use `q.displayln_info("CHECK: ...")` for intermediate test output — those lines are compared against reference `.log` files and make tests brittle.
- **Only one CHECK line**: Every test must end with `q.displayln_info("CHECK: finished successfully.")` as the final line.

## Code Style — Cython

### Directives (top of .pyx files)
```python
# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8
```

### Conventions
- `cimport` for C++ declarations, `import` for Python
- `snake_case` methods, `PascalCase` classes
- Section separators: `### -------------------------------------------------------------------`
- Error context: `fname = q.get_fname()` for error messages

## Key Directories

| Path | Contents |
|------|----------|
| `qlat-utils/qlat_utils/include/qlat-utils/` | C++ utility headers |
| `qlat/qlat/include/qlat/` | C++ core library headers |
| `qlat/qlat/lib/` | C++ source files |
| `qlat/qlat/*.pyx` | Cython bindings for the C++ core |
| `examples-py/` | Python test/example scripts — canonical reference; keep up to date |
| `examples-py-gpt/` | Python examples using the GPT/Grid interface — canonical reference; keep up to date |
| `examples-py-cps/` | Python examples using the CPS interface — canonical reference; keep up to date |
| `examples-cpp/` | C++ test/example programs — canonical reference; keep up to date |
| `examples-cpp-grid/` | C++ examples using the Grid interface — canonical reference; keep up to date |
| `scripts/` | Build and test shell scripts |
| `applications/` | **Old, possibly outdated code — do not reference** |

## Important Notes

- Do not commit changes to `.log` reference files unless intentionally updating test expectations
- MPI is required for most tests (`mpi4py`, `mpiexec`)
- CI (`.github/workflows/qlat.yml`) runs on push/PR to `master` and verifies no `CHECK:` lines changed
- Nix is available for reproducible builds (`nix-build` in `nixpkgs/`)
- Custom macros are prefixed `QLAT_` or `qacc_`
- Keep the `examples-*` directories up to date with API changes and run their tests routinely; they are the canonical, trustworthy reference for current usage
