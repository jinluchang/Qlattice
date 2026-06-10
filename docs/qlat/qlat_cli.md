# `qlat.__main__` — CLI Entry Point

Source: `qlat/qlat/__main__.py`

> **Note:** Update this document when updating the source file.

## Outline

- Usage message and subcommand dispatch.
- Delegates to `qlat_config.__main__` or `qlat.scripts.*` modules.

## Overview

The `__main__` module allows running qlat command-line tools via
`python3 -m qlat <subcommand>`.  Each subcommand imports and executes the
corresponding script module, which handles its own argument parsing and
MPI initialisation.

## Available Subcommands

| Subcommand | Script | Description |
|---|---|---|
| `qlat-config` | `qlat_config.__main__` | Print qlat build configuration |
| `eigen-system-checksum` | `scripts.eigen_system_checksum` | Verify eigen-system file checksums |
| `eigen-system-repartition` | `scripts.eigen_system_repartition` | Repartition eigen-system files |
| `fields-checksum` | `scripts.fields_checksum` | Compute checksums of field files |
| `fields-list` | `scripts.fields_list` | List field files in a directory |
| `fields-build-index` | `scripts.fields_build_index` | Build an index for field files |
| `fields-rewrite` | `scripts.fields_rewrite` | Rewrite field files (e.g. format change) |
| `fields-properly-truncate` | `scripts.fields_properly_truncate` | Truncate field files to valid data |
| `gauge-fix-coulomb` | `scripts.gauge_fix_coulomb` | Coulomb gauge-fix a gauge configuration |
| `topo-measure` | `scripts.topo_measure` | Measure topological charge |
| `topo-measure-wilson-flow` | `scripts.topo_measure_wilson_flow` | Measure topology via Wilson flow |
| `topo-measure-gpt` | `scripts.topo_measure_gpt` | Measure topology using GPT |

## Examples

### Print usage

```bash
python3 -m qlat
```

### Show qlat build configuration

```bash
python3 -m qlat qlat-config
```

### Checksum eigen-system files

```bash
python3 -m qlat eigen-system-checksum eigen-dir/
```

### List field files

```bash
python3 -m qlat fields-list fields-dir/
```

### Measure topological charge

```bash
python3 -m qlat topo-measure gauge-field-path
```
