# Coding doc structure

This document describes the standard procedure for adding documentation to Python
modules in `qlat-utils` and `qlat`.

## Overview

Each Python module (`.py` or `.pyx` file) should have:

1. A **module docstring** in the source file that references the documentation location.
2. A **markdown documentation file** in `docs/qlat-utils/` or `docs/qlat/`.

## Procedure

### Step 1: Identify the Module

Locate the source file to document. For example:
- `qlat-utils/qlat_utils/json.py`
- `qlat/qlat/field.pyx`

### Step 2: Add a Module Docstring to the Source File

Add a module docstring at the top of the source file (after any Cython
directives or shebang lines). The docstring must include:

1. The module name as a heading.
2. A one-line description of the module.
3. The path to the documentation file.
4. A note to update the documentation when updating the source.

**Template for `.py` files:**

```python
"""
Module ``package.module_name``
==============================

Short description of what this module provides.

Documentation: ``docs/package/module_name.md``

.. note:: Update the documentation when updating this source file.
"""
```

**Template for `.pyx` files (after the cython directive line):**

```python
# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``package.module_name``
==============================

Short description of what this module provides.

Documentation: ``docs/package/module_name.md``

.. note:: Update the documentation when updating this source file.
"""
```

If the module has existing attribution comments (e.g., license, original author),
move them into the docstring as a reference link at the bottom.

### Step 3: Create or Update the Documentation File

The documentation file lives in `docs/qlat-utils/` or `docs/qlat/` and is named
after the module (e.g., `qlat_json.md` for `qlat_utils/json.py`).

The document must begin with:

```markdown
# `package.module_name` — Short Description

Source: `package-path/module_name.py`

> **Note:** Update this document when updating the source file.
```

Then include:

1. **Outline** — A table of contents with anchor links.
2. **Overview** — What the module provides, key properties, usage summary.
3. **Detailed sections** — One section per major class, function group, or concept.
4. **Examples** — Practical usage examples at the end.

> **Important:** Code examples in `qlat` documentation must call
> `q.begin_with_mpi(size_node_list)` after imports to initialize the MPI
> environment before using any `qlat` functionality, and `q.end_with_mpi()`
> at the end to finalize properly. The `Geometry` constructor requires a
> `Coordinate` object (e.g., `q.Geometry(q.Coordinate([4, 4, 4, 8]))`).
> This is not required for `qlat-utils` modules that do not depend on MPI.

### Step 4: Add a Link in the Package Index

Add a `toctree` entry for the new documentation file in the appropriate package
index:

- `docs/qlat-utils.rst` for `qlat-utils` modules
- `docs/qlat.rst` for `qlat` modules

Add the entry in the `toctree` directive near the top of the file:

```rst
.. toctree::
   :maxdepth: 1

   qlat-utils/qlat_rng_state.md

   qlat-utils/qlat_json.md

   qlat-utils/qlat_new_module.md
```

### Step 5: Verify Examples

After adding or updating documentation that includes code examples, verify
that the examples actually work against the compiled module.

**Procedure:**

1. Create a fresh empty directory for testing:
    ```bash
    mkdir -p tmp/qlat-verify
    ```

2. Write a verification script that exercises all the code examples from the
   documentation, using ``assert`` statements to confirm correct behavior.
   The script should verify:
   - Round-trip fidelity for each supported type
   - Type preservation through nested structures
   - Edge cases (compact output, plain types, etc.)

 3. If the source code was modified, recompile with:
    ```bash
    bash scripts/qlat-all.sh
    ```

    If the build fails during the configuration step (e.g., ``meson`` error about
    missing dependencies or stale build state), perform a clean build first:
    ```bash
    bash scripts/qlat-clean-build.sh
    bash scripts/qlat-all.sh
    ```

4. Run the verification script from the fresh directory:
    ```bash
    cd tmp/qlat-verify && python3 verify_examples.py
    ```

**Example** — verifying ``qlat_utils.json``:

```bash
mkdir -p tmp/qlat-verify
# write tmp/qlat-verify/verify_json.py
python3 tmp/qlat-verify/verify_json.py
```

The verification script should print "All tests PASSED." on success. If any
assertion fails, fix the documentation or the source code, recompile, and
re-run until all examples pass.

### Step 6: Final Checks

Check that:
- The docstring renders correctly (e.g., `help(module)` or IDE hover).
- The markdown file is readable and all internal links work.
- The source path in the doc matches the actual file location.

## Naming Conventions

| Source File | Documentation File |
|---|---|
| `qlat_utils/json.py` | `docs/qlat-utils/qlat_json.md` |
| `qlat_utils/rng_state.pyx` | `docs/qlat-utils/qlat_rng_state.md` |
| `qlat/field.pyx` | `docs/qlat/qlat_field.md` |

The documentation filename uses the package prefix (e.g., `qlat_`) to avoid
collisions between packages.

