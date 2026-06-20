# `qlat.tempita` — Tempita Template Substitution

Source: `qlat/qlat/tempita.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Usage](#usage)
3. [Examples](#examples)

---

## Overview

The `qlat.tempita` module is a command-line utility for code generation
using [Tempita](https://github.com/cython/cython/blob/master/Cython/Tempita/_tempita.py)
template substitution, as bundled with Cython.

It reads a Tempita template file, performs `{{expr}}` substitution via
`Cython.Tempita.sub`, and writes the result to an output file. This is used
internally by the qlat build system to generate C++ / Cython source files from
parameterized templates.

The module is a standalone script (not an importable library). It takes exactly
two positional arguments — the input template path and the output file path.

This module is the qlat-package counterpart of `qlat_utils.tempita`. Both
are identical in functionality; they exist in separate packages so that each
package can invoke Tempita substitution without depending on the other.

---

## Usage

```bash
python -m qlat.tempita <template> <output>
```

| Argument | Description |
|---|---|
| `template` | Path to the Tempita template file |
| `output` | Path to write the generated output |

The template file is read as UTF-8, substituted, and the result is written as
UTF-8.

### Template Syntax

Tempita uses `{{expression}}` for inline substitution. Any Python expression
inside `{{ }}` is evaluated and its result is inserted into the output. For full
Tempita syntax (conditionals, loops, inheritance), see the Cython Tempita
documentation.

---

## Examples

### Basic Template Substitution

Given a template file `header.hpp.tempita`:

```
// Auto-generated — do not edit
#pragma once

constexpr int NDIM = {{ndim}};
constexpr double BETA = {{beta}};
```

Running:

```bash
python -m qlat.tempita header.hpp.tempita header.hpp
```

produces `header.hpp` with the `{{ndim}}` and `{{beta}}` placeholders replaced
by their substituted values.

### Build System Integration

In a Makefile or Meson build step, the module is typically invoked as:

```bash
python -m qlat.tempita template.pyx.tempita template.pyx
```

This pattern is used throughout the qlat codebase to generate Cython `.pyx` and
C++ `.hpp` files from parameterized templates.
