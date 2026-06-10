# `qlat_utils.cutils` — Low-Level Filesystem Utilities

Source: `qlat-utils/qlat_utils/cutils.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Path Manipulation](#path-manipulation)
3. [Directory Listing](#directory-listing)
4. [File and Directory Queries](#file-and-directory-queries)
5. [File and Directory Operations](#file-and-directory-operations)
6. [Cached Queries](#cached-queries)
7. [Info-Logged Operations](#info-logged-operations)
8. [MPI-Synchronized Operations](#mpi-synchronized-operations)
9. [Cache Management](#cache-management)
10. [Diagnostics](#diagnostics)
11. [Examples](#examples)

---

## Overview

`cutils` exposes the C++ filesystem layer of `qlat-utils` to Python. It provides:

- **Path manipulation** — `basename`, `dirname`, `all_dirname_vec`, `remove_trailing_slashes`
- **Directory listing** — `qls`, `qls_all` (local) and `qls_sync_node`, `qls_all_sync_node` (MPI-synchronized)
- **File/directory queries** — `does_file_exist`, `is_directory`, `is_regular_file`
- **File/directory operations** — `qmkdir`, `qmkdir_p`, `qrename`, `qremove`, `qremove_all`, `qtruncate`
- **Cached variants** — `_cache`-suffixed functions that memoize results for repeated queries
- **Info-logged variants** — `_info`-suffixed functions that log operations via `displayln_info`
- **MPI-synchronized variants** — `_sync_node`-suffixed functions that ensure all MPI ranks agree on filesystem state
- **Cache management** — `clear_all_caches`, `clear_mem_cache`, `get_all_caches_info`

Python access:

```python
import qlat_utils as q
```

---

## Path Manipulation

| Function | Description |
|---|---|
| `basename(fn)` | Return the filename portion of path `fn` (everything after the last `/`). |
| `dirname(fn)` | Return the directory portion of path `fn` (everything before the last `/`). |
| `all_dirname_vec(fn)` | Return a list of all ancestor directories of `fn`, from parent to root. |
| `remove_trailing_slashes(fn)` | Return `fn` with any trailing `/` characters removed. |

---

## Directory Listing

### `qls(path, is_sort=True) -> list[str]`

List the immediate contents of directory `path` (returned as full paths). If `is_sort` is `True`, results are sorted alphabetically.

### `qls_all(path, is_folder_before_files=False, is_sort=True) -> list[str]`

Recursively list all entries under `path` (returned as full paths). If `is_folder_before_files` is `True`, directories appear before their contents. If `is_sort` is `True`, entries are sorted within each directory level.

---

## File and Directory Queries

| Function | Return Type | Description |
|---|---|---|
| `does_file_exist(path)` | `bool` | `True` if `path` exists (file or directory). |
| `is_directory(path)` | `bool` | `True` if `path` is an existing directory. |
| `is_regular_file(path)` | `bool` | `True` if `path` is an existing regular file. |

---

## File and Directory Operations

| Function | Description |
|---|---|
| `qmkdir(path)` | Create a single directory. Fails if parent does not exist. |
| `qmkdir_p(path)` | Create directory and all missing parents (like `mkdir -p`). |
| `qrename(old_path, new_path)` | Rename/move `old_path` to `new_path`. |
| `qremove(path)` | Remove a single file or empty directory. |
| `qremove_all(path)` | Recursively remove `path` and all its contents (like `rm -rf`). |
| `qtruncate(path, offset=0)` | Truncate file at `path` to `offset` bytes (default 0 = empty). |

---

## Cached Queries

These functions memoize their results. Use `clear_is_directory_cache()` or `remove_entry_directory_cache(path)` to invalidate stale entries.

| Function | Description |
|---|---|
| `is_directory_cache(path)` | Cached version of `is_directory`. |
| `is_regular_file_cache(path)` | Cached version of `is_regular_file`. |
| `does_file_exist_cache(path)` | Cached version of `does_file_exist`. |
| `clear_is_directory_cache()` | Clear the entire directory-existence cache. |
| `remove_entry_directory_cache(path)` | Remove a single entry from the cache. |

---

## Info-Logged Operations

These are identical to their non-`_info` counterparts but log each operation via `displayln_info`, which is useful for debugging MPI runs.

| Function | Description |
|---|---|
| `qmkdir_info(path)` | `qmkdir` with logging. |
| `qmkdir_p_info(path)` | `qmkdir_p` with logging. |
| `qrename_info(old_path, new_path)` | `qrename` with logging. |
| `qremove_info(path)` | `qremove` with logging. |
| `qremove_all_info(path)` | `qremove_all` with logging. |

---

## MPI-Synchronized Operations

In an MPI environment, all ranks must agree on the filesystem state before proceeding. These `_sync_node` functions perform a barrier and broadcast so that every rank sees the same result.

| Function | Description |
|---|---|
| `qls_sync_node(path, is_sort=True)` | MPI-synchronized `qls`. |
| `qls_all_sync_node(path, is_folder_before_files=False, is_sort=True)` | MPI-synchronized `qls_all`. |
| `does_file_exist_sync_node(path)` | MPI-synchronized `does_file_exist`. |
| `is_directory_sync_node(path)` | MPI-synchronized `is_directory`. |
| `is_regular_file_sync_node(path)` | MPI-synchronized `is_regular_file`. |
| `qmkdir_sync_node(path)` | MPI-synchronized `qmkdir`. |
| `qmkdir_p_sync_node(path)` | MPI-synchronized `qmkdir_p`. |
| `qremove_sync_node(path)` | MPI-synchronized `qremove`. |
| `qremove_all_sync_node(path)` | MPI-synchronized `qremove_all`. |

---

## Cache Management

| Function | Description |
|---|---|
| `get_all_caches_info()` | Return a list of strings describing all internal caches and their sizes. |
| `clear_all_caches()` | Clear all internal caches (directory, mem, etc.). |
| `clear_mem_cache()` | Clear only the memory buffer cache. |

---

## Diagnostics

### `displayln_malloc_stats()`

Print glibc `malloc_stats()` output to stdout. Useful for tracking memory allocation in long-running jobs.

### `get_eigen_type() -> str`

Return the name of the Eigen backend type used by the C++ layer (e.g., `"Eigen"` or `"Grid"`).

---

## Examples

### Path Manipulation

```python
import qlat_utils as q

q.basename("/home/user/data/config.lime")   # "config.lime"
q.dirname("/home/user/data/config.lime")     # "/home/user/data"
q.all_dirname_vec("/a/b/c/d")                # ["/a/b/c", "/a/b", "/a", "/"]
q.remove_trailing_slashes("/home/user/")     # "/home/user"
```

### Listing and Querying

```python
import qlat_utils as q

entries = q.qls("/tmp/data")
print(entries)

if q.is_directory("/tmp/data/run-1"):
    print("exists")

if not q.does_file_exist("/tmp/data/output.bin"):
    print("not found")
```

### Creating and Removing Directories

```python
import qlat_utils as q

q.qmkdir_p("/tmp/data/run-1/gauge")
# ... do work ...
q.qremove_all("/tmp/data/run-1")
```

### MPI-Synchronized File Operations

```python
import qlat_utils as q

# All MPI ranks agree on whether the directory exists
if not q.is_directory_sync_node("/shared/output"):
    q.qmkdir_p_sync_node("/shared/output")
```

### Using Caches

```python
import qlat_utils as q

# Repeated queries hit the cache
for i in range(1000):
    if q.is_directory_cache("/tmp/data"):
        pass  # fast after first call

# Invalidate when filesystem changes
q.clear_is_directory_cache()
```

### Diagnostics

```python
import qlat_utils as q

info = q.get_all_caches_info()
for line in info:
    print(line)

q.displayln_malloc_stats()
print(q.get_eigen_type())
```
