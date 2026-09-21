# AGENTS.md — Guide for Coding Agents

Operating guide for coding agents working in this repository: project layout,
build and test commands, code style, documentation, lint and commit
conventions.

## Project Overview

Qlattice is a lattice QCD simulation library. The monorepo contains four
Python/C++ packages:

- `qlat-utils` — core utilities (C++17, Cython, Python): parallel-safe RNG
  (`RngState`), timers, coordinates, the QAR chunked file format, caches,
  `LatData`, data analysis (jackknife, fits, gnuplot wrapper).
- `qlat` — main lattice library (depends on qlat-utils): geometry, gauge and
  fermion fields, HMC, propagators, contractions (`qlat/auto_contractor/`),
  topology, Wilson flow, lattice I/O, job scripts.
- `qlat-grid` — Grid framework interface.
- `qlat-cps` — CPS framework interface.

`qlat` also installs the top-level `qlat_gpt` module (GPT interface), and it
re-exports everything from `qlat_utils`, so `import qlat as q` is the usual
entry point. The version is in `./VERSION`.

Primary languages: **C++17**, **Cython**, **Python**, **Bash**.

> **Note**: The `applications/` directory contains old code that may be outdated or no longer correct. Do NOT reference it as a source of truth or as an example to follow.
>
> **Note**: The `examples-*` directories (`examples-py/`, `examples-py-gpt/`, `examples-py-cps/`, `examples-cpp/`, `examples-cpp-grid/`) ARE the canonical reference for current API usage and are safe to reference. Keep them up to date when APIs change, and run their tests routinely to catch regressions.

## Build System

Build system: **Nix** which calls **Meson** to build the packages (via
`meson-python` for pip packaging). No CMake or npm. The same packages build
without Nix through `meson-python` (`pip install qlat-utils`,
`pip install qlat`), and the older `./build.sh` / `./scripts/*.sh` flow with
the `qcore/` compiler wrappers still works — `docs/contents/install.md`
documents both.

## Build and Testing

### Build qlat with nix

Use `nixpkgs/install-py-local-kernel-with-nix.py` to build qlat via nix. It creates a `./result-py-local` symlink to the nix store path. Build variants: `--variant cuda`, `--variant cudasupport`, `--variant cu`, `--variant clang`, `--variant pypi`, or `--all-variants` for every variant. See `--help` for details.

```bash
./nixpkgs/install-py-local-kernel-with-nix.py
```

A build takes a few minutes when it cannot reuse the nix cache, so run it in
the background while you do other work.

**Rebuild requirement**: nix snapshots the working tree when the build starts,
so any source edit invalidates it — edits made *while* a build runs are not
included. Before running a test, make sure `./result-py-local` exists and was
built from the current tree. To check a pure-Python module without rebuilding
(e.g. after a comment-only edit, or if you suspect a stale build):

```bash
diff qlat-utils/qlat_utils/data.py result-py-local/lib/python3.*/site-packages/qlat_utils/data.py
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

Each script copies sources into `./tmp/examples-*/` and runs the test there. After a run, check `./tmp/examples-*/<name>.py.p/` (Python/GPT/CPS) or `./tmp/examples-*/<name>/build/` (C++/Grid) for log files (`log.full.txt`, `log.txt`, `log.check.txt`).

Tests use **log-comparison**: each test prints `CHECK:` lines compared against reference `.log` files. See each `examples-*/Makefile` for the full test-running workflow.

**The runners never modify the original reference files.** They print a
reminder that the regenerated `.log`, `.log.json` and `.log.json.new` files
are only under `./tmp/examples-*/`; copy them into `examples-*/` yourself when
the change is intentional (see below).

### Add or update an example test

A test consists of the script plus two committed reference files. The `%.log`
rule in `examples-*/Makefile` runs the script with `--test`, which makes
`q.check_log_json()` compare every `q.json_results_append()` entry against
`<name>.log.json` (writing `<name>.log.json.new`), and extracts the
`CHECK: `/`INFO: `/`WARNING: ` lines into `<name>.log` (`.log.check.txt` is
just the `CHECK: ` lines and drives the pass/fail `diff`).

1. Write the test (see the Test Script Pattern below) and register it by
   adding `<name>.log` to the `tests = ...` list of `examples-*/Makefile`.
2. Run it once. With no reference yet, `q.check_log_json()` exits before the
   final `CHECK:` line is printed, so `log.txt` is empty and `make` reports
   `passed` vacuously: the first run only produces the reference content.
3. Install the generated reference and re-run, then install the log:

```bash
./nixpkgs/run-one-example-py.py <name>
cp tmp/examples-py/<name>.py.p/<name>.log.json.new examples-py/<name>.log.json
./nixpkgs/run-one-example-py.py <name>          # now runs to completion
cp tmp/examples-py/<name>.log examples-py/<name>.log
./nixpkgs/run-one-example-py.py <name>          # must print "passed"
```

After step 2, always check that
`tmp/examples-*/<name>.py.p/log.check.txt.new` is non-empty and ends with
`CHECK: finished successfully.` — an empty file means the script exited early
(a missing reference, an assertion, or a crash) and the `passed` line is
meaningless.

Reference `.log.json` files are compared relatively with `check_eps` (the
value passed to `q.check_log_json()`, e.g. `1e-10`), so last-bit differences
from a legitimate numerical change are tolerated, but a changed or added entry
name is a mismatch. Regenerate a reference only when the change is intended —
never edit reference logs to make a failing test pass.

### Run all tests

**REQUIREMENT**: All tests MUST be run via the `nixpkgs/build-many-qlat-pkgs.py` script below. Do NOT loop through individual tests manually — the `tests` package set builds the `qlat-tests` package for the default version and name, which runs the whole test suite (Python, GPT, CPS, C++ and Grid examples) and fails when any `log.check.txt` differs:

```bash
./nixpkgs/build-many-qlat-pkgs.py --group tests
```

The script sets up the nix and nom caches (falling back to `<repo>/tmp/` when the user directories are not writable) and writes the test package to `$HOME/qlat-build/nix/tests/result` (`<repo>/tmp/qlat-build/nix/tests/result` when `$HOME` is not writable). `--group` may be repeated and combined with `-j`/`--cores`, and `./nixpkgs/build-many-qlat-pkgs.py --list-groups` shows every package set (`tests`, `core`, `all`, `cuda*`, `small`, ...).

**DO NOT** use shell loops like `for test in ... ; do ./nixpkgs/run-one-example-py.py $test ; done` — this is incorrect and bypasses the proper test orchestration.

### Run a new program with the nix-built environment

Always run new programs under `<project_root>/tmp/` (or a sub-directory) to keep test artifacts out of the source tree.

Source `setenv-qlat.sh` from the nix result to configure PATH, PYTHONPATH, LD_LIBRARY_PATH, etc. (`bind-gpu-qlat.sh` and `mpiexec` are then on `PATH`):

```bash
source result-py-local/bin/setenv-qlat.sh
mkdir -p tmp
cd tmp
mpiexec -n 2 --oversubscribe --bind-to none bash bind-gpu-qlat.sh python3 -m mpi4py ./your-script.py
```

See the `%.log` rule in `examples-py/Makefile` for the full invocation pattern (mpi options, CHECK-line extraction, etc.).

### Fast iteration on pure-Python changes

For a quick edit/run loop, shadow the installed package instead of rebuilding.
Site-packages entries are symlinks into the nix store, so copy with `-L`:

```bash
source result-py-local/bin/setenv-qlat.sh
mkdir -p tmp/pylib
cp -rL "$(python3 -c 'import qlat_utils, os; print(os.path.dirname(qlat_utils.__file__))')" tmp/pylib/qlat_utils
chmod -R u+w tmp/pylib/qlat_utils
cp qlat-utils/qlat_utils/data.py tmp/pylib/qlat_utils/data.py
cd tmp
PYTHONPATH=$PWD/pylib mpiexec -n 2 --oversubscribe --bind-to none python3 -m mpi4py ./script.py
```

The compiled extensions still come from the store build, so this validates
Python logic only. A rebuild plus the official runner is always required for
the final verification. `result-py-local` is just a symlink, so it can also be
repointed at another build temporarily (e.g. to check an earlier commit) as
long as it is restored afterwards.

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

Inside the packages themselves a module keeps its Python level qlat imports in
a module local `class q` block (as in `qlat_utils/data.py` and
`qlat_utils/jackknife_utils.py`) and refers to them as `q.<name>`, so that a
module does not re-export what it merely imports; non qlat imports
(`import numpy as np`, `import sys`, ...) stay at module level, and the same
holds for `cimport` in a `.pyx`. The exception is a module whose whole purpose
is to re-export — `qlat_utils/c.py`, the `qlat_grid`/`qlat_cps` `c.py`
aggregators and their `init.py`/`prop.py` facades — which keep their
`from .x import *` / re-export list, and imports deliberately placed after the
definitions to break an import cycle, which keep their position. When moving
such an import, keep the package `__init__` and the `__all__` lists in mind:
they must not lose a name.

### Conventions
- Shebang: `#!/usr/bin/env python3` for executable scripts
- `snake_case` for functions/variables, `PascalCase` for classes
- f-strings for formatting
- No type hints in example scripts (but welcomed in library code)
- No blank lines inside function bodies: use an indented `#` separator line
  instead. `lint.py` enforces this (see below).

### MPI conventions
- Start with `q.begin_with_mpi()` / `q.begin_with_gpt()` /
  `q.begin_with_grid()` and finish with the matching `q.end_with_mpi()` /
  `q.end_with_gpt()` / `q.end_with_grid()`.
- Node identity: `q.get_id_node()` / `q.get_num_node()`. For collective
  operations use the communicator from `q.get_comm()`, which satisfies
  `comm.rank == q.get_id_node()` (the qlat node order differs from
  `MPI_COMM_WORLD` for some Grid layouts).
- Every node must call collective operations the same number of times; guard
  node-specific work with `if q.get_id_node() == 0:` only for things like
  printing and single-node file writes.
- Split work over nodes with `q.get_distributed_range(total, id_node,
  num_node)` (contiguous, node-ordered, covers the whole range) and use
  `q.get_collective_comm(tag)` for new collective MPI code.

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
- **Comparisons**: `q.check_log_json(__file__, check_eps=...)` compares the
  recorded floats relatively; when a test also asserts something itself, use a
  tolerance (e.g. a max-relative-difference helper) rather than
  `np.array_equal` unless exact equality is really part of the contract.

### Data analysis
Jackknife and error analysis has a dedicated, documented API: use
`q.g_mk_jk()` / `q.g_jk_avg_err()` / `q.g_jk_size()` with the shared settings
in `q.default_g_jk_kwargs` (or the `q.JkKwargs(...)` context manager) instead
of hand-rolling jackknife code. Collective MPI variants are
`q.g_mk_jk(..., is_sync_node=True)` (same input on every node) and
`q.g_mk_jk_distributed()` (input split between the nodes). See
`docs/qlat-utils/qlat_jackknife_utils.md` and
`docs/contents/how_to_analysis.md`.

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

## Lint and Formatting

`lint.py` is the repository formatter; run it on the files you touch before
committing (the tree is expected to be lint-clean):

```bash
python3 lint.py path/to/file.py        # rewrite in place
python3 lint.py qlat-utils/qlat_utils  # a whole directory
python3 lint.py --check path/to/file.py   # exit 1 if changes are needed
python3 lint.py --diff path/to/file.py    # show the proposed changes
```

It formats `.py` with `ruff`, `.pyx` with a regex-based pass and C/C++
`.c/.cc/.cpp/.cxx/.h/.hpp` with `clang-format`, and it replaces empty lines
inside function/method bodies with indented comment lines (`#` for
Python/Cython, `//` for C/C++). Module-level blank lines and blank lines
between functions are left alone.

## Documentation

Each module's docstring names its documentation file and asks you to update
it:

```
Documentation: ``docs/qlat-utils/qlat_data.md``
.. note:: Update the documentation when updating this source file.
```

85 modules currently follow this convention. The mapping is `docs/<package>/qlat_<module>.md`
(e.g. `qlat-utils/qlat_utils/data.py` → `docs/qlat-utils/qlat_data.md`,
`qlat/qlat/hmc.pyx` → `docs/qlat/qlat_hmc.md`); when in doubt, follow the path
in the docstring. `docs/contents/` holds the prose guides (installation,
environment variables, command line arguments, analysis workflow, ...) and
`docs/*.rst` / `docs/*.md` cover the larger subsystems (jackknife,
auto-contractor, `qlat_gpt`, `qlat-scripts-v1`, ...).

The Sphinx site is built from `docs/` and deployed to
<https://jinluchang.github.io/Qlattice> by CI, so:
- document new public functions, modules and CLI options in the matching file;
- keep the guides consistent with the code when behavior changes;
- use the same reStructuredText-ish style as the surrounding text.

## Key Directories

| Path | Contents |
|------|----------|
| `qlat-utils/qlat_utils/` | Python modules and Cython sources of `qlat_utils` |
| `qlat-utils/qlat_utils/include/qlat-utils/` | C++ utility headers |
| `qlat/qlat/` | Python/Cython sources of `qlat` |
| `qlat/qlat/include/qlat/` | C++ core library headers |
| `qlat/qlat/lib/` | C++ source files |
| `qlat/qlat_scripts/` | Job and analysis scripts (e.g. `v1/`) |
| `qlat/auto_contractor/` | Automatic contraction code generation |
| `qlat/qlat_gpt.py` | GPT interface, installed as a top-level module |
| `qlat-grid/`, `qlat-cps/` | Grid and CPS interface packages |
| `examples-py/` | Python test/example scripts — canonical reference; keep up to date |
| `examples-py-gpt/` | Python examples using the GPT/Grid interface — canonical reference; keep up to date |
| `examples-py-cps/` | Python examples using the CPS interface — canonical reference; keep up to date |
| `examples-cpp/` | C++ test/example programs — canonical reference; keep up to date |
| `examples-cpp-grid/` | C++ examples using the Grid interface — canonical reference; keep up to date |
| `docs/` | Sphinx documentation: guides in `contents/`, per-module docs in `qlat-utils/`, `qlat/`, ... |
| `nixpkgs/` | Nix expressions and the build/test scripts described above |
| `scripts/` | Per-machine build scripts and dependency installers (non-Nix flow) |
| `qcore/` | Compiler wrappers and `setenv.sh` used by `scripts/` |
| `tmp/` | Scratch space for builds and test runs (git-ignored) |
| `distfiles/`, `projects/`, `cmd-*` | Local, git-ignored downloads, project setups and personal shortcuts — do not commit or rely on them |
| `applications/` | **Old, possibly outdated code — do not reference** |

## Important Notes

- Do not commit changes to `.log` / `.log.json` reference files unless intentionally updating test expectations (see "Add or update an example test")
- MPI is required for most tests (`mpi4py`, `mpiexec`)
- CI (`.github/workflows/qlat.yml`) runs on push/PR to `master`: it builds `nixpkgs/q-pkgs.nix` attributes (`pkgs-std-ucxless.qlat-env`/`qlat-tests`, then the GPT/CPS and jhub-env variants), so any change must survive a full nix build and the whole test suite; it also publishes the docs and the PyPI packages
- Nix is available for reproducible builds (`nix-build` in `nixpkgs/`); the result links (`result*`) and `tmp/` are git-ignored
- Custom macros are prefixed `QLAT_` or `qacc_`
- Keep the `examples-*` directories up to date with API changes and run their tests routinely; they are the canonical, trustworthy reference for current usage
- `release.py` performs the release: build the tarballs, tag, push, create the GitHub release, upload to PyPI and bump `VERSION`

## Commits

- Commit in small, topic-grouped commits with a short imperative subject and a
  detailed body explaining *what* changed and *why* (match the existing `git log`
  style); reference the affected modules and mention any intentional reference
  log updates.
- Keep unrelated changes out of the same commit; a large change that spans
  several concerns is easier to review as one commit per concern.
- The working tree must be clean and the tests must pass before committing.
