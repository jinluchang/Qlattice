# `qlat_utils.c` — Central C Extension Loader and API Aggregator

Source: `qlat-utils/qlat_utils/c.py`

> **Note:** Update this document when updating the source file.

## Outline

- Native shared-library loading (`libqlat-utils.so` / `.dylib`)
- Re-export of all public symbols from submodules
- Coordinate operations (integer and floating-point)
- RNG state management
- Matrix types: `WilsonMatrix`, `SpinMatrix`, `ColorMatrix`
- Element types for typed field containers
- MPI / node utilities
- Timer and profiling system
- File-system helpers with caching and sync variants
- QAR archive I/O
- Lattice data containers (`LatData`)

## Overview

`qlat_utils.c` is the single entry point that **loads the compiled C++ shared
library** into the Python process and then **re-exports every public symbol**
from the internal submodules (`timer`, `cutils`, `types`, `coordinate`,
`lat_data`, `rng_state`, `qar`, `mat`).  Importing this module is equivalent
to importing the entire low-level C-backed API of `qlat-utils`.

Before loading the library the module temporarily enables `RTLD_GLOBAL` so
that dependent Cython extensions can resolve symbols from the same shared
object.

## Detailed Sections

### Shared-Library Loading

```python
# Automatically locates and loads libqlat-utils.so (or .dylib on macOS)
# Uses RTLD_GLOBAL so that Cython extensions link against the same object.
```

The path is resolved relative to the package directory under `lib/`.

### Coordinate Operations

Functions for mapping between lattice coordinates and linear indices, and for
modular / signed-modular arithmetic used in momentum-space and periodic
boundary conditions:

| Function | Description |
|----------|-------------|
| `Coordinate`, `CoordinateD` | Integer and double coordinate types |
| `mod_coordinate(c, L)` | `c mod L` (standard modular) |
| `smod_coordinate(c, L)` | Signed modular: result in `[-L/2, L/2)` |
| `smod_sym_coordinate(c, L)` | Symmetric signed modular |
| `middle_mod_coordinate(c, L)` | Midpoint modular |
| `coordinate_from_index(index, size)` | Map linear index to coordinate tuple |
| `index_from_coordinate(coordinate, size)` | Map coordinate tuple to linear index |

The `_d` variants operate on double-precision coordinates.

### RNG State

`RngState` is a reproducible pseudo-random number generator with a seed
hierarchy.  `get_data_sig` computes a deterministic signature of a data
buffer, and `random_permute` produces a random permutation.

### Matrix Types and Operations

Three fixed-size matrix types used in lattice QCD spin-color algebra:

- **`WilsonMatrix`** — 12×12 complex matrix (spin ⊗ color)
- **`SpinMatrix`** — 4×4 complex matrix (spin space only)
- **`ColorMatrix`** — 3×3 complex matrix (color space only)

Key operations follow a naming convention `mat_<op>_<type1>_<type2>`:

| Prefix | Operation |
|--------|-----------|
| `mat_tr_` | Trace |
| `mat_mul_` | Multiplication |
| `mat_add_` | Addition |

Scalar multiplication uses `mat_mul_a_*`.  Gamma matrices are accessed via
`get_gamma_matrix`.  `wilson_matrix_g5_herm` computes the g5-hermitian
conjugate, and `mat_epsilon_contraction_wm_wm_wm` performs the epsilon
contraction over three Wilson matrices.

### Element Types

`ElemType*` constants describe the scalar type stored in field containers
(e.g. `ElemTypeComplexD`, `ElemTypeWilsonMatrix`, `ElemTypeRealF`, etc.).
They are used by the higher-level `Field` and `SelectedField` classes.

### MPI / Node Utilities

| Function | Description |
|----------|-------------|
| `get_id_node()` | MPI rank of the current process |
| `get_num_node()` | Total number of MPI processes |
| `sync_node()` | Global MPI barrier |

### Timer and Profiling

A hierarchical, scope-based timer system:

| Function / Class | Description |
|------------------|-------------|
| `Timer` | Context-manager timer (enabled by default) |
| `TimerNone` | No-op timer for disabled timers |
| `timer(name)` | Decorator / context manager |
| `timer_verbose(name)` | Timer that prints on destruction |
| `timer_flops(name, flops)` | Timer with FLOP count for MFLOP/s reporting |
| `get_time()` | Current wall-clock time |
| `get_total_time()` | Elapsed time since start |
| `get_remaining_time()` | Time left before the configured limit |
| `timer_display()` | Print accumulated timer statistics |
| `timer_reset()` | Reset all timer counters |

`displayln` and `displayln_info` are the standard output functions; the latter
prepends MPI rank information.

### File-System Helpers

Wrapper functions that add logging (`_info` suffix), MPI-synchronous variants
(`_sync_node` suffix), and an `is_directory` cache layer:

- `qls`, `qls_all` — directory listing
- `qmkdir`, `qmkdir_p`, `qrename`, `qremove`, `qremove_all`
- `does_file_exist`, `is_directory`, `is_regular_file`
- Cache-aware: `is_directory_cache`, `clear_is_directory_cache`

### QAR Archive I/O

QAR is a simple tar-like archive format used by Qlattice:

- `qar_create`, `qar_extract`, `qar_build_index`
- `QarFile`, `QFile` — file-handle wrappers
- `qcat`, `qcat_bytes`, `qtouch`, `qappend` — convenience I/O
- `compute_crc32` — CRC-32 checksum

### Lattice Data (`LatData`)

`LatData` (and variants `LatDataRealF`, `LatDataInt`, `LatDataLong`) is a
ranked array with named dimensions, commonly used for per-site or per-momentum
data.  `mk_lat_data` creates a new instance; `load_lat_data` reads one from
disk.

## Examples

### Load the library and inspect available symbols

```python
import qlat_utils as q

# The module is loaded; all sub-symbols are available
print("Number of MPI nodes:", q.get_num_node())
print("Node ID:", q.get_id_node())
```

### Use coordinate helpers

```python
import qlat_utils as q

size = (4, 4, 4, 8)
coord = q.coordinate_from_index(13, size)
print("Index 13 ->", coord)           # (1, 0, 0, 1)
idx = q.index_from_coordinate(coord, size)
print("Back to index:", idx)          # 13

print("smod(-1, 4) =", q.smod_coordinate(-1, 4))  # 3
```

### Matrix operations

```python
import qlat_utils as q

gamma = q.get_gamma_matrix(0)        # gamma_0 SpinMatrix
wm = q.WilsonMatrix()                # zero 12x12 matrix
tr = q.mat_tr_wm(wm)                # trace
wm2 = q.mat_mul_wm_wm(wm, wm)       # matrix multiply
wm_herm = q.wilson_matrix_g5_herm(wm)  # g5-hermitian conjugate
```

### Timer usage

```python
import qlat_utils as q

with q.Timer("my-section"):
    # timed block
    pass

@q.timer_verbose("important-function")
def expensive():
    pass
```

### QAR archive

```python
import qlat_utils as q

# Create an archive from a directory
q.qar_create("output.qar", "data_dir", False)

# Extract it
q.qar_extract("output.qar", "extracted_dir", False)
```
