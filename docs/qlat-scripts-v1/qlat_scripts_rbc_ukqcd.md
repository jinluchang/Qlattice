# `qlat_scripts.v1.rbc_ukqcd` — RBC/UKQCD Fermion Actions, Eigensolvers, and Inverters

Source: `qlat/qlat_scripts/v1/rbc_ukqcd.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Parameter Retrieval](#parameter-retrieval)
3. [Fermion Action and Preconditioners](#fermion-action-and-preconditioners)
4. [Eigensolvers](#eigensolvers)
5. [Eigensystem I/O](#eigensystem-io)
6. [Inverter Construction](#inverter-construction)
7. [Examples](#examples)

---

## Overview

This module provides the interface between Qlattice and GPT for RBC/UKQCD lattice QCD ensembles. It handles:

- Fermion parameter lookup with automatic fallback to lower accuracy levels
- Domain-wall fermion action setup (Mobius and ZMobius)
- Even-odd preconditioning (eo1, eo2, kappa-ne)
- Lanczos eigensolver and Chebyshev-accelerated deflated solver construction
- Eigensystem save/load (uncompressed and compressed formats)
- Mixed-precision CG inverter construction with deflation, split-grid, and MADWF support
- A native Qlattice inverter path for specific ensembles (24D, 32D)

All functions use `get_param(job_tag, ...)` from `rbc_ukqcd_params` to read ensemble-specific configuration.

## Terminology

| Term | Meaning |
|------|---------|
| `job_tag` | Ensemble identifier, e.g. `"24D"`, `"64I"`, `"test-4nt8"` |
| `inv_type` | Quark flavor index: `0` = light, `1` = strange |
| `inv_acc` | Inversion accuracy: `0` = sloppy, `1` = medium, `2` = exact |
| `gf` | Qlattice gauge field (`q.GaugeField`) |
| `gpt_gf` | GPT gauge field (from `qg.gpt_from_qlat`) |
| `eig` | Eigensystem for deflation (`EigSystemGPT` or `EigSystemCompressedGPT`) |
| `pc_ne` | Even-odd preconditioned operator |
| `parity` | Checkerboard parity (typically `g.odd`) |

---

## Parameter Retrieval

### `get_param_fermion(job_tag, inv_type, inv_acc)`

Returns fermion action parameters for the given ensemble and quark flavor. Falls back to lower `inv_acc` values if the requested accuracy level is not configured.

**Returns:** `dict` or `None`.

### `get_ls_from_fermion_params(fermion_params)`

Extracts the fifth-dimensional extent `Ls` from fermion parameters. For ZMobius actions, `Ls` is inferred from the length of the `"omega"` list; otherwise it is read from the `"Ls"` key.

**Returns:** `int`.

### `get_param_lanc(job_tag, inv_type, inv_acc=0)`

Returns Lanczos parameters with fallback to lower `inv_acc`.

**Returns:** `dict` or `None`.

### `get_param_clanc(job_tag, inv_type, inv_acc=0)`

Returns Chebyshev-Lanczos (coarse-grid) parameters with fallback to lower `inv_acc`.

**Returns:** `dict` or `None`.

### `get_param_cg_mp_maxiter(job_tag, inv_type, inv_acc)`

Returns the CG mixed-precision max iteration count. Reads from `cg_params-{inv_type}-{inv_acc}` if configured; otherwise uses defaults: `50` for test ensembles, `500` for `inv_acc >= 2`, `200` for light quarks, `300` for strange.

**Returns:** `int`.

---

## Fermion Action and Preconditioners

### `mk_pc_parity(job_tag, inv_type, inv_acc)`

Returns the checkerboard parity used for even-odd preconditioning. Currently always returns `g.odd`.

### `mk_pc_ne(job_tag, inv_type, inv_acc, *, eig=None, parity=None)`

Constructs the even-odd preconditioned operator. Selects the preconditioner type based on the ensemble and fermion parameters:

- `eo2_kappa_ne` — for ZMobius actions (when `"omega"` is in fermion params and no eigenvectors are provided); prints a warning as this does not support split-CG deflation
- `eo1_ne` — for ensembles `"64I"` and `"64I-pq"`
- `eo2_ne` — default for all other ensembles

**Returns:** GPT preconditioner callable.

### `mk_quark_matrix(job_tag, gpt_gf, inv_type, inv_acc)`

Creates the quark matrix (fermion operator) for the given gauge field. Uses `g.qcd.fermion.zmobius` if the fermion parameters contain `"omega"`, otherwise `g.qcd.fermion.mobius`.

**Returns:** GPT quark matrix object.

---

## Eigensolvers

### `mk_eig(job_tag, gf, inv_type, inv_acc=0, *, parity=None, pc_ne=None)`

Computes eigenvectors of the preconditioned fermion operator using implicitly restarted Lanczos (IRL) with Chebyshev polynomial acceleration.

**Algorithm:**
1. Constructs the quark matrix and preconditioner
2. Runs power iteration to estimate the largest eigenvalue
3. Runs IRL with Chebyshev polynomial filtering
4. Returns an `EigSystemGPT` object containing eigenvectors and eigenvalues

**Requires:** `get_param(job_tag, "lanc_params", inv_type, inv_acc)`.

**Returns:** `qg.EigSystemGPT`.

### `mk_ceig(job_tag, gf, inv_type, inv_acc=0, *, parity=None, pc_ne=None)`

Computes compressed (deflated) eigenvectors via coarse-grid Chebyshev-Lanczos. Falls back to `mk_eig` if no `clanc_params` are configured.

**Algorithm:**
1. Computes fine-grid eigenvectors via IRL (same as `mk_eig`)
2. Constructs a coarse-grid subspace from the leading `nbasis` eigenvectors
3. Builds the coarse operator and runs IRL on it
4. Smooths coarse eigenvectors with a CG solver
5. Returns an `EigSystemCompressedGPT` object

**Requires:**
- `get_param(job_tag, "lanc_params", inv_type, inv_acc)`
- `get_param(job_tag, "clanc_params", inv_type, inv_acc)`

**Returns:** `qg.EigSystemCompressedGPT`.

### `get_smoothed_evals(basis, cevec, gf, job_tag, inv_type, inv_acc=0, *, parity=None, pc_ne=None)`

Recomputes smoothed eigenvalues for an existing compressed eigensystem. Useful when the smoother parameters change but the basis and coarse eigenvectors remain valid.

**Returns:** `list[float]` — smoothed eigenvalues.

---

## Eigensystem I/O

### `save_eig(path, eig, job_tag, inv_type=0, inv_acc=0)`

Saves an eigensystem to disk. For compressed eigensystems (`EigSystemCompressedGPT`), uses `clanc_params["save_params"]` to control single-precision storage and MPI I/O. For uncompressed eigensystems (`EigSystemGPT`), saves directly.

### `load_eig_lazy(path, job_tag, inv_type=0, inv_acc=0)`

Lazily loads an eigensystem from disk. Returns `None` if the path does not exist or required files are missing. Otherwise returns a `@q.lazy_call`-wrapped callable that loads the eigensystem on first invocation.

Detects two formats:
- **Compressed** — expects `metadata.txt`, `eigen-values.txt`, and either `00/000000000000.compressed` or `00.zip`
- **Uncompressed** — expects `global`, `index`, and `index.crc32`

**Returns:** `callable` (returns `EigSystemCompressedGPT` or `EigSystemGPT`) or `None`.

---

## Inverter Construction

### `mk_gpt_inverter(gf, job_tag, inv_type, inv_acc, *, gt=None, mpi_split=None, n_grouped=None, eig=None, eps=1e-8, parity=None, pc_ne=None, qtimer=True)`

Constructs a GPT-based mixed-precision inverter with full feature support.

**Features:**
- **Deflation** — uses eigenvectors (compressed or uncompressed) to accelerate convergence
- **Split-grid** — distributes inversions across sub-communicators via `mpi_split`
- **Grouped inversions** — batches multiple sources via `n_grouped`
- **MADWF** — mixed-action domain-wall fermion deflation when `Ls` differs between sloppy and exact solvers
- **Mixed precision** — single/double precision defect-correcting CG
- **Gauge transformation** — wraps inverter with optional gauge transform

**Parameters:**
- `gt` — gauge transformation for twisted-boundary-condition inverter
- `mpi_split` — MPI sub-communicator layout, e.g. `[1, 1, 1, 4]`
- `n_grouped` — number of sources solved simultaneously
- `eig` — eigensystem for deflation
- `eps` — solver residual tolerance
- `qtimer` — timer object or `True`/`False` for auto/no timer

**Returns:** `qg.InverterGPT` or `q.InverterGaugeTransform`.

### `mk_qlat_inverter(gf, job_tag, inv_type, inv_acc, *, gt=None)`

Constructs a native Qlattice domain-wall inverter. Only supported for ensembles `"24D"` and `"32D"`.

**Returns:** `q.InverterDomainWall` or `q.InverterGaugeTransform`.

### `mk_inverter(*args, **kwargs)`

Convenience alias for `mk_gpt_inverter`.

### `get_inv(gf, job_tag, inv_type, inv_acc, *, gt=None, mpi_split=None, n_grouped=None, eig=None, eps=1e-8, pc_ne=None, qtimer=True)`

Cached inverter factory. Checks `q.cache_inv` for an existing inverter matching the given parameters before constructing a new one. Avoids redundant inverter builds when the same gauge field and configuration are reused.

**Returns:** `qg.InverterGPT` or `q.InverterGaugeTransform`.

---

## Examples

### Basic Inverter Setup

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])
import qlat_scripts.v1 as qs

job_tag = "24D"
inv_type = 0   # light quark
inv_acc = 2    # exact solve

# gf and gt assumed available from prior setup
inv = qs.rbc_ukqcd.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt)

q.end_with_mpi()
```

### Eigensolver with Deflation

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])
import qlat_scripts.v1 as qs

job_tag = "24D"
inv_type = 0

# Compute eigenvectors (fine-grid Lanczos)
eig = qs.rbc_ukqcd.mk_eig(job_tag, gf, inv_type)

# Build inverter with deflation
inv = qs.rbc_ukqcd.get_inv(gf, job_tag, inv_type, inv_acc=0, eig=eig)

# Save eigenvectors for reuse
qs.rbc_ukqcd.save_eig("eig-cache", eig, job_tag, inv_type)

q.end_with_mpi()
```

### Compressed Eigenvectors

```python
import qlat as q
q.begin_with_mpi([[1, 1, 1, 4]])
import qlat_scripts.v1 as qs

job_tag = "24D"
inv_type = 0

# Compute compressed eigenvectors (coarse-grid deflation)
ceig = qs.rbc_ukqcd.mk_ceig(job_tag, gf, inv_type)

# Load lazily from disk
load_fn = qs.rbc_ukqcd.load_eig_lazy("eig-cache", job_tag, inv_type)
if load_fn is not None:
    eig = load_fn()

q.end_with_mpi()
```

---

## Key Parameters

These parameters are read via `get_param(job_tag, ...)`:

| Parameter | Description |
|-----------|-------------|
| `total_site` | Lattice dimensions, e.g. `[24, 24, 24, 64]` |
| `fermion_params` | Fermion action parameters per `(inv_type, inv_acc)` |
| `cg_params-{inv_type}-{inv_acc}` | CG solver parameters (`maxiter`, `maxcycle`, `pv_maxiter`) |
| `lanc_params` | Lanczos parameters (`fermion_params`, `pit_params`, `cheby_params`, `irl_params`) |
| `clanc_params` | Chebyshev-Lanczos parameters (`nbasis`, `block`, `cheby_params`, `irl_params`, `smoother_params`, `save_params`) |

---

## Notes

- All `@q.timer_verbose` decorated functions provide automatic timing and memory instrumentation.
- The `@q.lazy_call` decorator in `load_eig_lazy` ensures eigensystem loading happens at most once per callable.
- MADWF (mixed-action domain-wall fermion) is automatically enabled when the sloppy and exact solvers use different `Ls` values. It is disabled for strange quarks (`inv_type == 1`) to avoid using eigenvectors for the exact solve.
- The `mk_qlat_inverter` path is a fallback for ensembles without GPT support and only handles `"24D"` and `"32D"`.
- Split-grid inversions (`mpi_split`) distribute work across MPI sub-communicators, improving throughput for large ensembles.
