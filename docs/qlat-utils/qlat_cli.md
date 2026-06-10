# `qlat_utils.__main__` — Command-Line Interface

Source: `qlat-utils/qlat_utils/__main__.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Usage](#usage)
3. [Subcommands](#subcommands)
   - [qlat-utils-config](#qlat-utils-config)
   - [crc32](#crc32)
   - [lat-io-glimpse](#lat-io-glimpse)
   - [lat-io-diff](#lat-io-diff)
   - [pickle-glimpse](#pickle-glimpse)
   - [qar-glimpse](#qar-glimpse)
   - [qar](#qar)
4. [Examples](#examples)

---

## Overview

The `qlat_utils.__main__` module provides the command-line interface (CLI) entry
point for the `qlat_utils` package. It is invoked via:

```bash
python3 -m qlat_utils <subcommand> [args...]
```

The module dispatches to one of several subcommands for inspecting, comparing, and
managing lattice QCD data files. Each subcommand is implemented in a separate script
under `qlat_utils/scripts/`.

---

## Usage

```
python3 -m qlat_utils qlat-utils-config ...
python3 -m qlat_utils crc32 ...
python3 -m qlat_utils lat-io-glimpse ...
python3 -m qlat_utils lat-io-diff ...
python3 -m qlat_utils pickle-glimpse ...
python3 -m qlat_utils qar-glimpse ...
python3 -m qlat_utils qar ...
```

Running without arguments prints the usage message and exits.

---

## Subcommands

### `qlat-utils-config`

```bash
python3 -m qlat_utils qlat-utils-config [args...]
```

Delegates to the `qlat_utils_config.__main__` module for querying build
configuration (compiler flags, library paths, etc.).

### `crc32`

```bash
python3 -m qlat_utils crc32 <file1> [file2] ...
```

Compute and print the CRC32 checksum of one or more files.

**Output format:** One line per file: `<hash> '<filename>'`

The hash is computed in 64 KB chunks using `zlib.crc32` and printed as an 8-digit
hexadecimal string.

### `lat-io-glimpse`

```bash
python3 -m qlat_utils lat-io-glimpse <file1> [file2] ...
```

Load and display the contents of one or more `LatData` files (`.lat` format).

Each file is loaded via `q.LatData().load(fn)` and its contents are printed to
stdout with a header line showing the filename.

### `lat-io-diff`

```bash
python3 -m qlat_utils lat-io-diff <file1> <file2>
```

Compare two `LatData` files and print the norm of their difference.

**Output (three lines):**
1. `qnorm(ld1 - ld2)` — squared norm of the difference.
2. `qnorm(ld1)` — squared norm of the first file.
3. `qnorm(ld2)` — squared norm of the second file.

Exactly two filenames must be provided.

### `pickle-glimpse`

```bash
python3 -m qlat_utils pickle-glimpse <file1> [file2] ...
```

Load and pretty-print the contents of one or more pickle files.

Each file is loaded via `q.load_pickle_obj(fn)` and displayed using
`pprint.pformat` with a header line showing the filename.

### `qar-glimpse`

```bash
python3 -m qlat_utils qar-glimpse <path1.qar> [path2.qar] ...
```

List the contents of one or more QAR (Qlattice Archive) files.

Recursively enumerates entries in `.qar` archives, printing each entry with a
sequential index, source archive path, internal index, and the relative path of the
entry. Nested `.qar` files are traversed automatically.

### `qar`

```bash
python3 -m qlat_utils qar <action> [args...]
```

Full QAR archive management tool. Supports actions such as creating, extracting,
listing, and verifying QAR archives. This is a more feature-rich version of
`qar-glimpse` with additional operations.

---

## Examples

### Computing CRC32 Checksums

```bash
python3 -m qlat_utils crc32 data/config_001.lat data/config_002.lat
```

Output:
```
a1b2c3d4 'data/config_001.lat'
e5f6a7b8 'data/config_002.lat'
```

### Inspecting a Lattice Data File

```bash
python3 -m qlat_utils lat-io-glimpse data/propagator.lat
```

### Comparing Two Lattice Data Files

```bash
python3 -m qlat_utils lat-io-diff data/prop_v1.lat data/prop_v2.lat
```

Output:
```
1.23456e-10
4.56789e+02
4.56789e+02
```

### Inspecting a Pickle File

```bash
python3 -m qlat_utils pickle-glimpse results/analysis.pkl
```

### Listing QAR Archive Contents

```bash
python3 -m qlat_utils qar-glimpse data/gauge_config.qar
```

Output:
```
       0 'data/gauge_config.qar'        0 'header.txt'
       1 'data/gauge_config.qar'        1 'field.lat'
       2 'data/gauge_config.qar'        2 'metadata.json'
```
