# AGENTS.md — Guide for Coding Agents

## Project Overview

Qlattice is a lattice QCD simulation library. The monorepo contains four Python/C++ packages:
- `qlat-utils` — core utilities (C++17, Cython, Python)
- `qlat` — main lattice library (depends on qlat-utils)
- `qlat-grid` — Grid framework interface
- `qlat-cps` — CPS framework interface

Primary languages: **C++17**, **Cython**, **Python**, **Bash**.

## Build System

Build system: **Meson** (via `meson-python` for pip packaging). No CMake or npm.

### Install via custom scripts
```bash
./scripts/setenv.default.sh
./scripts/qcore.sh
./scripts/qlat-utils.sh
./scripts/qlat.sh
./scripts/qlat-grid.sh
./scripts/qlat-cps.sh
```

### Full build orchestration
```bash
./nixpkgs/install-py-local-kernel-with-nix.sh
```

## Testing

Tests live in `examples-py/` (Python) and `examples-cpp/` (C++). The testing framework is **log-comparison**: each test prints `CHECK:` lines, which are compared against committed reference `.log` files.

### Run all tests
```bash
./nixpkgs/build-many-qlat-pkgs-core.sh -j 4 --cores 15
```

### Run a single Python test
```bash
cd examples-py
make utils.log                        # runs utils.py, compares CHECK lines
```
Test execution pattern (from Makefile):
```bash
mpiexec -n 2 --oversubscribe --bind-to none python3 -m mpi4py ./utils.py \
  --test -qmp-geom 1 1 1 2 --mpi 1.1.1.2 --mpi_split 1 1 1 1 \
  --mpi 1.1.2 --mpi_split 1 1 1
```

### Run a single C++ test
```bash
cd examples-cpp
make simple-1.run                     # builds and runs simple-1
```

### Verify test output (CI check)
```bash
! git diff | grep '^[+-]CHECK: '     # fails if CHECK lines changed
```

### Update reference logs after intentional changes
```bash
cd examples-py && make run           # regenerate logs
git add *.log                        # commit new reference logs
```

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
- Use project macros: `qassert(cond)`, `Qassert(cond)`, `qerr(msg)`
- Warnings: `warn(msg)` or `displayln_info(ssprintf("WARNING: ..."))`
- Test markers: `displayln_info("CHECK: ...")` for test verification

### Key Macros/Decorators
- `qacc` — accelerator decorator (CPU/GPU portability)
- `TIMER("name")` — scope-based timer
- `API` — DLL export macro

## Code Style — Python

### Imports
```python
import sys
import qlat as q
import qlat_utils as qu   # if needed
```
Standard library first, then `qlat` (aliased as `q`).

### Conventions
- Shebang: `#!/usr/bin/env python3` for executable scripts
- `snake_case` for functions/variables, `PascalCase` for classes
- f-strings for formatting
- No type hints in example scripts (but welcomed in library code)

### Test Script Pattern
```python
q.begin_with_mpi(size_node_list)
q.json_results_append(f"test description")
# ... test logic ...
q.json_results_append(f"{result}")
```
Every test must end with a `CHECK: finished successfully.` line (via `q.displayln_info`).

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

## Formatting & Linting

- **C++**: Run `clang-format -i` (uses `.clang-format` config)
- **Python/Cython**: No enforced formatter. Follow existing style (4-space indent, no trailing whitespace)
- No flake8, ruff, black, or pylint configs exist — do not add them

## Key Directories

| Path | Contents |
|------|----------|
| `qlat-utils/qlat_utils/include/qlat-utils/` | C++ utility headers |
| `qlat/qlat/include/qlat/` | C++ core library headers |
| `qlat/qlat/lib/` | C++ source files |
| `qlat/cqlat/` | Cython bindings |
| `examples-py/` | Python test/example scripts |
| `examples-cpp/` | C++ test/example programs |
| `scripts/` | Build and test shell scripts |

## Important Notes

- Do not commit changes to `.log` reference files unless intentionally updating test expectations
- MPI is required for most tests (`mpi4py`, `mpiexec`)
- CI (`.github/workflows/qlat.yml`) runs on push/PR to `master` and verifies no `CHECK:` lines changed
- Nix is available for reproducible builds (`nix-build` in `nixpkgs/`)
- Custom macros are prefixed `QLAT_` or `qacc_`
