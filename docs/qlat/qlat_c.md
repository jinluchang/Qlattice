# `qlat.c` — Central C Extension Loader and API Aggregator

Source: `qlat/qlat/c.py.in`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Shared-Library Loading](#shared-library-loading)
3. [Exported API Categories](#exported-api-categories)
4. [Re-exported Submodules](#re-exported-submodules)
5. [Examples](#examples)

---

## Overview

`qlat.c` is the single entry point that **loads the compiled C++ shared
library** (`libqlat.so` / `libqlat.dylib`) into the Python process and then
**re-exports every public symbol** from the internal Cython and Python
submodules.  Importing this module is equivalent to importing the entire
low-level C-backed API of `qlat`.

The module is generated from a Mako template (`.py.in`) that expands
type-specific `Field{Type}`, `SelectedField{Type}`, and `SelectedPoints{Type}`
class names for every lattice data type.

Before loading the library the module temporarily enables `RTLD_GLOBAL` so
that dependent Cython extensions can resolve symbols from the same shared
object.

---

## Shared-Library Loading

```python
# Automatically locates and loads libqlat.so (or libqlat.dylib on macOS)
# Uses RTLD_GLOBAL so that Cython extensions link against the same object.
```

The path is resolved relative to the package directory under `lib/`.

---

## Exported API Categories

The module's `__all__` list groups symbols by functional area:

| Category | Key Symbols |
|----------|-------------|
| **MPI / Node** | `begin`, `end`, `get_size_node`, `get_coor_node`, `bcast_*`, `glb_sum_*` |
| **Geometry** | `Geometry`, `geo_resize`, `geo_eo` |
| **Fields** | `Field{Type}`, `SelectedField{Type}`, `SelectedPoints{Type}`, `Field`, `FieldBase` |
| **Field Selection** | `PointsSelection`, `FieldSelection`, `mk_xg_field`, `get_psel_tslice` |
| **Shuffle** | `SelectedShufflePlan`, `split_fields`, `merge_fields` |
| **Lock / Quit** | `obtain_lock`, `release_lock`, `qquit`, `check_time_limit`, `check_stop` |
| **Field Operations** | `field_expanded`, `refresh_expanded`, `CommPlan`, `mk_phase_field`, `FastFourierTransform`, `field_shift`, `shuffle_field` |
| **Gauge Field** | `GaugeField`, `GaugeTransform`, `gf_avg_plaq`, `gf_avg_link_trace`, `gf_plaq_field` |
| **Propagators** | `Prop`, `SelProp`, `PselProp`, `SpinProp`, `FermionField4d`, `set_point_src`, `free_invert` |
| **Smearing** | `gf_ape_smear`, `gf_hyp_smear`, `prop_smear` |
| **Topology** | `gf_topology_field`, `gf_topology`, `gf_topology_terms` |
| **Wilson Flow** | `gf_wilson_flow_force`, `gf_wilson_flow_step`, `gf_wilson_flow`, `gf_energy_density` |
| **Gauge Action** | `GaugeAction` |
| **HMC** | `GaugeMomentum`, `set_rand_gauge_momentum`, `gf_evolve`, `run_hmc_pure_gauge` |
| **Muon Line** | `calc_muon_line_m`, `get_muon_line_m`, `compute_save_muonline_interpolation` |
| **HLBL Contract** | `mk_m_z_field_tag`, `contract_four_pair_labels`, `mk_local_current_from_props` |
| **Fields I/O** | `ShuffledFieldsWriter`, `ShuffledFieldsReader`, `open_fields`, `list_fields`, `fields_build_index` |

---

## Re-exported Submodules

After loading the shared library, the module imports and re-exports symbols
from these internal submodules:

| Submodule | Domain |
|-----------|--------|
| `qlat_utils.c` | All qlat-utils low-level symbols |
| `cqlat` | Cython bindings for qlat C++ core |
| `utils_io` | File I/O utilities |
| `mpi` | MPI communication |
| `geometry` | Lattice geometry |
| `field_base` | Field base classes |
| `field_types` | Typed field classes |
| `field_selection` | Field/point selection |
| `selected_field_types` | Selected field classes |
| `selected_points_types` | Selected points classes |
| `field_utils` | Field utility operations |
| `qcd` | QCD-specific operations |
| `propagator` | Propagator types and operations |
| `smear` | Smearing operations |
| `hmc` | Hybrid Monte Carlo |
| `gauge_action` | Gauge action |
| `topology` | Topological charge |
| `wilson_flow` | Wilson flow |
| `muon_line` | Muon line integrals |
| `hlbl_contract` | HLBL contractions |
| `fields_io` | Shuffled fields I/O |

---

## Examples

### Load the library and inspect available symbols

```python
import qlat as q

q.begin_with_mpi([])

print("Number of MPI nodes:", q.get_num_node())
print("Node ID:", q.get_id_node())

q.end_with_mpi()
```

### Create a geometry and field

```python
import qlat as q

q.begin_with_mpi([])

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))
gf = q.GaugeField(geo)
gf.set_unit()

print("Average plaq:", q.gf_avg_plaq(gf))

q.end_with_mpi()
```

### Matrix operations

```python
import qlat_utils as q

gamma = q.get_gamma_matrix(0)
wm = q.WilsonMatrix()
tr = q.mat_tr_wm(wm)
wm2 = q.mat_mul_wm_wm(wm, wm)
```
