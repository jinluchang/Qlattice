# How To: Writing Analysis Code with Qlattice

## Project Structure

```
project/
├── run.py                  # Data generation (measurement) code
├── analysis/
│   └── my_analysis.py      # Analysis scripts
├── results/                # Data output directory
│   └── {job_tag}/
│       └── {measurement}/
│           └── traj-{traj}/
│               └── *.lat   # LatData files
└── how_to_analysis.md      # This file
```

## 1. Setup

### Imports

```python
import glob
import os

import numpy as np

import qlat as q
import qlat_scripts.v1 as qs
```

### Load Path Configuration

Always set `load_path_list` at module level before any data loading:

```python
qs.load_path_list[:] = [
    "results",
    "../results",
]
```

### MPI Lifecycle

Every qlattice script must start/end with MPI:

```python
if __name__ == "__main__":
    q.begin_with_mpi()
    # ... analysis code ...
    q.end_with_mpi()
```

## 2. Discovering Available Data

### Finding Trajectories

Use `qs.get_load_path()` to resolve paths across multiple search directories:

```python
@q.cache_call(maxsize=128)
@q.timer
def get_traj_list(job_tag):
    path = qs.get_load_path(f"{job_tag}/my_measurement")
    traj_list = []
    if path is None:
        return traj_list
    pattern = os.path.join(path, "traj-*")
    for d in glob.glob(pattern):
        basename = os.path.basename(d)
        if os.path.isdir(d) and basename.startswith("traj-"):
            try:
                traj = int(basename[5:])
                traj_list.append(traj)
            except ValueError:
                pass
    traj_list.sort()
    return traj_list
```

**Key points:**
- `qs.get_load_path(fn)` returns the first matching path or `None`
- Always check for `None` before using the path
- Sort trajectory lists for reproducibility

### Finding Files for a Trajectory

```python
fn = f"{job_tag}/my_measurement/traj-{traj}/data.lat"
path = qs.get_load_path(fn)
assert path is not None, f"Cannot find {fn}"
```

## 3. Loading Data

### LatData Format

LatData files (`.lat`) store multi-dimensional arrays with metadata:

```python
ld = q.LatData()
ld.load(path)

# Access metadata
info = ld.info()
# info[i] = [name, size, [labels...], ...]

# Access data as numpy array
data = ld[:]  # shape: (dim0, dim1, ...)
```

### Example: Loading Correlator Data

```python
@q.timer
def load_data(job_tag, traj):
    fn = f"{job_tag}/measurement/traj-{traj}/data.lat"
    path = qs.get_load_path(fn)
    assert path is not None
    ld = q.LatData()
    ld.load(path)
    info = ld.info()
    
    # Extract expression names from first dimension
    expr_names = info[0][2]  # list of str
    
    # Extract sizes from remaining dimensions
    sizes = [info[i][1] for i in range(1, len(info))]
    
    # Get data array
    data = ld[:].copy()  # always copy to avoid issues
    
    return expr_names, data, sizes
```

## 4. Data Processing

### Efficient Binning with numpy

Use `np.add.at` for unbuffered accumulation (unlike `+=` which is buffered):

```python
@q.timer
def bin_by_coordinate(data, coord_arr):
    """
    Bin complex data by coordinate values.
    
    Parameters
    ----------
    data : np.ndarray, shape (n_expr, *spatial_dims), dtype complex
    coord_arr : np.ndarray, shape (*spatial_dims), dtype int
    
    Returns
    -------
    coord_values : np.ndarray, sorted unique coordinate values
    binned : np.ndarray, shape (n_expr, n_bins), dtype complex
    """
    coord_flat = coord_arr.ravel()
    coord_values = np.unique(coord_flat)
    max_coord = int(coord_values.max()) + 1
    n_expr = data.shape[0]
    
    # Sum per bin for each expression
    data_flat = data.reshape(n_expr, -1)
    sums = np.zeros((n_expr, max_coord), dtype=complex)
    for e in range(n_expr):
        np.add.at(sums[e], coord_flat, data_flat[e])
    
    # Keep only valid bins
    valid = np.bincount(coord_flat, minlength=max_coord) > 0
    binned = sums[:, valid]
    coord_values = np.where(valid)[0]
    
    return coord_values, binned
```

**Key points:**
- Use **sum** not mean for jackknife analysis (jackknife needs raw sums)
- `np.add.at` handles repeated indices correctly
- Filter out empty bins with `np.bincount`

### Building Coordinate Arrays

```python
# For a grid with dimensions [size_t, size_z, size_y, size_x]
coords = [np.arange(s) for s in sizes]
grids = np.meshgrid(*coords, indexing="ij")
x_sqr = sum(g**2 for g in grids)  # x^2 + y^2 + z^2 + t^2
```

## 5. Jackknife Error Estimation

### Using `q.g_mk_jk`

The qlattice jackknife utilities support two modes: randomized jackknife (`"rjk"`, the default) and Super-Jackknife (`"super"`). The default `"rjk"` mode is recommended for most use cases:

```python
@q.timer
def collect_and_jackknife(job_tag, traj_list):
    """
    Load data for all trajectories and apply jackknife.
    
    Returns
    -------
    metadata : tuple
        Expression names, coordinate values, etc.
    jk_arr : np.ndarray, shape (1 + n_rand_sample, n_expr, n_bins)
        Index 0 is the mean; indices 1..n are jackknife replicates.
    """
    all_data = []
    jk_idx_list = []
    metadata = None
    
    for traj in traj_list:
        expr_names, data, sizes = load_data(job_tag, traj)
        if metadata is None:
            metadata = (expr_names, sizes)
        
        coord_values, binned = bin_by_coordinate(data, coord_arr)
        all_data.append(binned)
        jk_idx_list.append((job_tag, traj))  # IMPORTANT: use (job_tag, traj) tuples
    
    samples = np.array(all_data)  # shape: (n_traj, n_expr, n_bins)
    
    # Apply Super-Jackknife
    jk_arr = q.g_mk_jk(samples, jk_idx_list)
    
    return metadata, jk_arr
```

### Extracting Results

**Do not extract avg/err until right before printing/using:**

```python
jk_arr = collect_and_jackknife(job_tag, traj_list)

# Extract avg and err only when needed
avg = q.g_jk_avg(jk_arr)
err = q.g_jk_err(jk_arr)

# Print results
print(f"Mean: {avg[0, 0]:.6e} +/- {err[0, 0]:.6e}")
```

**Key points:**
- `jk_idx_list` must be `[(job_tag, traj), ...]` tuples, not integers
- `jk_arr[0]` is the mean; `jk_arr[1:]` are replicates
- Keep `jk_arr` around for further analysis (plotting, combining datasets)
- `q.g_mk_jk` accepts numpy arrays directly (no `.tolist()` needed)
- Use default `n_rand_sample` unless you have a reason to change it

### Available jk Utilities

```python
q.g_mk_jk(data_list, jk_idx_list)   # Create jackknife array
q.g_jk_avg(jk_arr)                   # Extract mean (index 0)
q.g_jk_err(jk_arr)                   # Extract error estimate
q.g_jk_avg_err(jk_arr)               # Extract both at once
```

## 6. Output Formatting

### Printing Tables

```python
def format_results(expr_names, coord_values, avg, err, max_display=10):
    """Print results: value +/- error."""
    # Header
    print(f"{'coord':>8s}", end="")
    for name in expr_names:
        print(f"  {name:>30s}", end="")
    print()
    
    # Data rows
    for k, x in enumerate(coord_values):
        if x > max_display:
            break
        print(f"{x:8d}", end="")
        for e in range(len(expr_names)):
            v = avg[e, k]
            s = err[e, k]
            print(f"  ({v.real:12.6e} +/- {s.real:12.6e})", end="")
        print()
```

## 7. Complete Template

```python
#!/usr/bin/env python3

import glob
import os

import numpy as np

import qlat as q
import qlat_scripts.v1 as qs

qs.load_path_list[:] = [
    "results",
    "../results",
]


@q.cache_call(maxsize=128)
@q.timer
def get_traj_list(job_tag):
    """Discover available trajectories for a measurement."""
    path = qs.get_load_path(f"{job_tag}/my_measurement")
    traj_list = []
    if path is None:
        return traj_list
    pattern = os.path.join(path, "traj-*")
    for d in glob.glob(pattern):
        basename = os.path.basename(d)
        if os.path.isdir(d) and basename.startswith("traj-"):
            try:
                traj = int(basename[5:])
                traj_list.append(traj)
            except ValueError:
                pass
    traj_list.sort()
    return traj_list


@q.timer
def load_data(job_tag, traj):
    """Load LatData file and return (names, data, metadata)."""
    fn = f"{job_tag}/my_measurement/traj-{traj}/data.lat"
    path = qs.get_load_path(fn)
    assert path is not None, f"Cannot find {fn}"
    ld = q.LatData()
    ld.load(path)
    info = ld.info()
    names = info[0][2]
    data = ld[:].copy()
    return names, data


@q.timer
def process_data(data):
    """Process raw data into analysis format."""
    # Your processing here
    # Example: return np.real(data)
    return data


@q.timer
def collect_data(job_tag, traj_list):
    """Load all trajectories, process, and apply jackknife."""
    all_data = []
    jk_idx_list = []
    names = None
    
    for traj in traj_list:
        n, raw = load_data(job_tag, traj)
        if names is None:
            names = n
        processed = process_data(raw)
        all_data.append(processed)
        jk_idx_list.append((job_tag, traj))
    
    samples = np.array(all_data)
    jk_arr = q.g_mk_jk(samples, jk_idx_list)
    return names, jk_arr


def format_results(names, avg, err):
    """Print formatted results table."""
    for i, name in enumerate(names):
        print(f"{name}: {avg[i].real:.6e} +/- {err[i].real:.6e}")


if __name__ == "__main__":
    q.begin_with_mpi()
    
    job_tag = "48I"
    traj_list = get_traj_list(job_tag)
    print(f"{job_tag}: {len(traj_list)} trajectories")
    
    if len(traj_list) > 0:
        names, jk_arr = collect_data(job_tag, traj_list)
        avg = q.g_jk_avg(jk_arr)
        err = q.g_jk_err(jk_arr)
        format_results(names, avg, err)
    
    q.end_with_mpi()
```

## 8. Common Pitfalls

| Mistake | Fix |
|---------|-----|
| Forgetting `q.begin_with_mpi()` / `q.end_with_mpi()` | Always wrap main in MPI lifecycle |
| Using `+=` for binning | Use `np.add.at` for unbuffered accumulation |
| Converting to mean before jackknife | Keep sums; jackknife needs raw data |
| Using integer indices for `jk_idx_list` | Use `(job_tag, traj)` tuples |
| Extracting avg/err inside loop | Extract once at the end from `jk_arr` |
| Calling `.tolist()` before `g_mk_jk` | Pass numpy arrays directly |
| Not checking `get_load_path` return | Always check for `None` |
| Forgetting to copy LatData | Use `ld[:].copy()` to avoid reference issues |

## 9. Performance Tips

- Use `@q.timer` on all functions to profile
- Use `@q.cache_call(maxsize=N)` for expensive lookups
- Prefer vectorized numpy operations over Python loops
- Use `np.add.at` for binning (not `np.bincount` + indexing)
- Keep `jk_arr` in memory for multiple analyses (don't recompute)
