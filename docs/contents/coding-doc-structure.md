# Coding doc structure

This document describes the standard procedure for adding documentation to Python
modules in `qlat-utils`, `qlat`, `qlat-scripts-v1`, and `auto-contractor`.

## Overview

Each Python module (`.py`, `.pyx`, or `.py.in`/`.pyx.in` template) should have:

1. A **module docstring** in the source file that references the documentation location.
2. A **markdown documentation file** in the appropriate `docs/` subdirectory:

| Package | Source | Docs Directory |
|---------|--------|---------------|
| qlat-utils | `qlat-utils/qlat_utils/` | `docs/qlat-utils/` |
| qlat | `qlat/qlat/` | `docs/qlat/` |
| qlat-scripts-v1 | `qlat/qlat_scripts/v1/` | `docs/qlat-scripts-v1/` |
| auto-contractor | `qlat/auto_contractor/` | `docs/auto-contractor/` |

## Procedure

### Step 1: Identify the Module

Locate the source file to document. For example:
- `qlat-utils/qlat_utils/json.py`
- `qlat/qlat/field.pyx`
- `qlat/qlat/field_types.pyx.in` (template that generates `.pyx`)
- `qlat/qlat/c.py.in` (template that generates `.py`)
- `qlat/qlat_gpt.py`
- `qlat/qlat_scripts/v1/gen_data.py`
- `qlat/auto_contractor/operators.py`

Note: `qlat/qlat/*.py.in` and `qlat/qlat/*.pyx.in` are template files that
generate the actual `.py`/`.pyx` modules. Document the generated module, not
the template itself.

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

**Template for `.py.in` files (after the Mako `{{py: ...}}` block):**

```python
{{py:
# Mako template definitions ...
}}

"""
Module ``package.module_name``
==============================

Short description of what this module provides.

Documentation: ``docs/package/module_name.md``

.. note:: Update the documentation when updating this source file.
"""
```

**Template for `.pyx.in` files (after the cython directive and Mako block):**

```python
# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

{{py:
# Mako template definitions ...
}}

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

The documentation file lives in the appropriate `docs/` subdirectory and is named
after the module. Examples:
- `docs/qlat-utils/qlat_json.md` for `qlat_utils/json.py`
- `docs/qlat/qlat_field.md` for `qlat/field.pyx`
- `docs/qlat_gpt.md` for `qlat/qlat_gpt.py` (top-level package module)
- `docs/qlat-scripts-v1/qlat_scripts_gen_data.md` for `qlat_scripts/v1/gen_data.py`
- `docs/auto-contractor/auto_contractor_operators.md` for `auto_contractor/operators.py`

Top-level package modules (files directly under `qlat/` such as `qlat_gpt.py`)
place their documentation at the `docs/` root (e.g., `docs/qlat_gpt.md`) rather
than under a package subdirectory. Their toctree entry goes in `docs/index.rst`
rather than `docs/qlat.rst`.

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
>
> Code that is intended to run on a single process (e.g., file I/O, printing,
> or computations that do not involve MPI collectives) should be wrapped in
> ``if q.get_id_node() == 0:`` to avoid duplicate output or conflicting writes
> when running under MPI with multiple processes.
>
> For `auto-contractor` modules that use GPT, code examples must call
> `qg.begin_with_gpt()` after imports and `qg.end_with_gpt()` at the end.

### Step 4: Add a Link in the Package Index

Add a `toctree` entry for the new documentation file in the appropriate package
index:

- `docs/qlat-utils.rst` for `qlat-utils` modules
- `docs/qlat.rst` for `qlat` modules (files under `qlat/qlat/`)
- `docs/index.rst` for top-level `qlat` package modules (e.g., `qlat/qlat_gpt.py`)
- `docs/qlat-scripts-v1.rst` for `qlat-scripts-v1` modules
- `docs/auto-contractor.rst` for `auto-contractor` modules

Add the entry in the `toctree` directive near the top of the file:

```rst
.. toctree::
   :maxdepth: 1

   qlat-utils/qlat_rng_state.md

   qlat-utils/qlat_json.md

   qlat-utils/qlat_new_module.md
```

For `qlat-scripts-v1` or `auto-contractor`, the pattern is the same but with
the corresponding directory name:

```rst
.. toctree::
   :maxdepth: 1

   qlat-scripts-v1/qlat_scripts_gen_data.md

   auto-contractor/auto_contractor_operators.md
```

### Step 5: Verify Examples

After adding or updating documentation that includes code examples, verify
that the examples actually work against the compiled module.

**Procedure:**

1. Build the environment with nix:
    ```bash
    cd nixpkgs && name='' ./install-py-local-kernel-with-nix.sh
    ```

2. Source the environment:
    ```bash
    source result-py-local/bin/setenv-qlat.sh
    ```

3. Write a verification script in the appropriate examples directory.

   - For modules that **do not** require GPT/Grid, write to `examples-py/`:
     `examples-py/docs-<module_name>.py`

   - For modules that **require** GPT/Grid, write to `examples-py-gpt/`:
     `examples-py-gpt/docs-<module_name>.py`

   The script should exercise all the code examples from the documentation,
   using ``assert`` statements to confirm correct behavior. The script should
   verify:
   - Round-trip fidelity for each supported type
   - Type preservation through nested structures
   - Edge cases (compact output, plain types, etc.)

4. Run the verification script using the standard test runner:
    ```bash
    # For non-GPT scripts:
    ./nixpkgs/run-one-example-py.py docs-<module_name>

    # For GPT scripts:
    ./nixpkgs/run-one-example-py-gpt.py docs-<module_name>
    ```

**Example** — verifying ``qlat_utils.json``:

```bash
# write examples-py/docs-qlat_utils_json.py
./nixpkgs/run-one-example-py.py docs-qlat_utils_json
```

The verification script should print "CHECK: finished successfully." at the
end (this is the standard test completion marker used by the test runner). If
any assertion fails, fix the documentation or the source code, rebuild, and
re-run until all examples pass.

5. Add the new test to the `tests` (or `tests_gpt`) list in the corresponding
   Makefile so it is included in the full test suite:

   - `examples-py/Makefile` — add `docs-<module_name>.log \` to the `tests` list.
   - `examples-py-gpt/Makefile` — add `docs-<module_name>.log \` to the `tests_gpt` list.

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
| `qlat/field_types.pyx.in` | `docs/qlat/qlat_field_types.md` |
| `qlat/selected_field_types.pyx.in` | `docs/qlat/qlat_selected_field_types.md` |
| `qlat/selected_points_types.pyx.in` | `docs/qlat/qlat_selected_points_types.md` |
| `qlat/c.py.in` | `docs/qlat/qlat_c.md` |
| `qlat/qlat_gpt.py` | `docs/qlat_gpt.md` |
| `qlat_scripts/v1/gen_data.py` | `docs/qlat-scripts-v1/qlat_scripts_gen_data.md` |
| `qlat_scripts/v1/rbc_ukqcd.py` | `docs/qlat-scripts-v1/qlat_scripts_rbc_ukqcd.md` |
| `auto_contractor/operators.py` | `docs/auto-contractor/auto_contractor_operators.md` |
| `auto_contractor/auto_contract_compilation.py` | `docs/auto-contractor/auto_contractor_auto_contract_compilation.md` |

Template files (`.py.in`, `.pyx.in`) use Mako templating to generate the
actual `.py`/`.pyx` modules. Document the **template file** as the source
(e.g., `Source: qlat/qlat/field_types.pyx.in`) since that is what developers
edit. The existing docs for `field_types`, `selected_field_types`, and
`selected_points_types` already follow this convention.

The documentation filename uses the package prefix (e.g., `qlat_`) to avoid
collisions between packages. For `qlat-scripts-v1` modules, the source lives in
`qlat_scripts/v1/` and documentation uses the `qlat_scripts_` prefix. For
`auto-contractor` modules, the documentation uses the `auto_contractor_` prefix.

