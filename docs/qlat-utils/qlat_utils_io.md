# `qlat_utils.utils_io` — File I/O, Hashing, and Caching Utilities

Source: `qlat-utils/qlat_utils/utils_io.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Directory Helpers](#directory-helpers)
3. [Byte I/O](#byte-io)
4. [JSON I/O](#json-io)
5. [Pickle I/O](#pickle-io)
6. [Hashing](#hashing)
7. [Caching Decorators](#caching-decorators)
8. [Display Method](#display-method)
9. [Examples](#examples)

---

## Overview

The `qlat_utils.utils_io` module provides file I/O utilities with MPI-aware
synchronisation, JSON and pickle serialisation wrappers, SHA-256 hashing for
cache keys, and decorators for transparent result caching — both in-memory
(LRU) and on-disk.

All I/O operations are exposed as free functions under the `qlat_utils` (`q`)
namespace via the internal `class q` trick, so they are accessed as
`q.save_json_obj(...)`, `q.load_pickle_obj(...)`, etc.

---

## Directory Helpers

### `qmkdirs(path)`

Create directories at `path` (including parents), then flush the directory
cache for the entry.

```python
q.qmkdirs("results/run_001/data")
```

### `qmkdirs_info(path)`

Same as `qmkdirs` but prints a log message on node 0 before creating directories.

### `mk_dirs(path)` / `mk_dirs_info(path)`

Lower-level alternatives identical in behaviour to `qmkdirs` / `qmkdirs_info`.

### `mk_file_dirs(fn)`

Create the parent directory of a file path `fn` (so that `fn` can be written).
Does nothing if the file is in the current directory.

```python
q.mk_file_dirs("output/data.pickle")
```

### `mk_file_dirs_info(path)`

Same as `mk_file_dirs` but prints a log message on node 0.

---

## Byte I/O

### `save_bytes(obj_str, path, *, is_sync_node=True)`

Save a `bytes` or `str` object to `path`. When `is_sync_node=True` (default),
only node 0 writes; all other nodes skip the write. Automatically calls
`mk_file_dirs_info(path)`.

```python
q.save_bytes(b"hello world", "data.bin")
```

### `load_bytes(path, default_value=None, *, is_sync_node=True)`

Read a file as `bytes`. When `is_sync_node=True` (the default), a single node
reads the file and broadcasts to all other nodes via the QMP/Grid runtime.
When `is_sync_node=False`, every node reads independently.

Returns `default_value` if the file is absent.

```python
data = q.load_bytes("data.bin")
data = q.load_bytes("missing.bin", default_value=b"fallback")
```

---

## JSON I/O

### `save_json_obj(obj, path, *, indent=2, is_sync_node=True)`

Serialize `obj` to JSON and save to `path`. Uses `json_dumps` internally, which
handles NumPy and extended numeric types (see `qlat_json.md`).

Set `indent=None` for compact (single-line) output.

```python
q.save_json_obj({"mass": 0.5, "iterations": 1000}, "params.json")
```

### `load_json_obj(path, default_value=None, *, is_sync_node=True)`

Load a JSON file and deserialize with `json_loads`, restoring extended types.

```python
params = q.load_json_obj("params.json")
params = q.load_json_obj("missing.json", default_value={})
```

---

## Pickle I/O

### `save_pickle_obj(obj, path, *, is_sync_node=True)`

Pickle `obj` and save to `path` as raw bytes.

```python
q.save_pickle_obj(my_array, "data.pickle")
```

### `load_pickle_obj(path, default_value=None, *, is_sync_node=True)`

Load and unpickle from `path`. Returns `default_value` if the file does not
exist.

```python
arr = q.load_pickle_obj("data.pickle")
arr = q.load_pickle_obj("missing.pickle", default_value=[])
```

---

## Hashing

### `hash_sha256(s)`

Compute the SHA-256 hex digest of `s`. Supports:

| Input Type | Behaviour |
|---|---|
| `str` | Encoded as UTF-8, then hashed |
| `bytes` | Hashed directly |
| `tuple` | Hashes `"tuple:"` prefix + hash of each element |
| `list` | Hashes `"list:"` prefix + hash of each element |
| `np.ndarray` | Hashes `"np.ndarray:"` + `repr(tolist())` |
| custom object | Delegates to `s.hash_sha256()` if the attribute exists |
| other | Falls back to `repr(s)` |

```python
key1 = q.hash_sha256("hello")
key2 = q.hash_sha256((3.14, "lattice", 42))
key3 = q.hash_sha256(np.array([1.0, 2.0, 3.0]))
```

---

## Caching Decorators

### `pickle_cache_call(func, path, *, is_sync_node=True)`

Call `func()` once and cache the result on disk at `path` as a pickle file.
On subsequent calls the cached result is returned.

If `is_sync_node=True` (default), all nodes see the same cached value.

```python
def expensive_computation():
    # heavy work
    return result

result = q.pickle_cache_call(expensive_computation, "cache/result.pickle")
```

### `pickle_cache(path, is_sync_node=True)` (decorator)

Decorator version of `pickle_cache_call`. The path is a directory; the function
arguments (pickled and hashed) determine the cache filename within that directory.

```python
@q.pickle_cache("cache/my_func")
def my_func(x, y):
    return x + y

result = my_func(3, 4)       # computes and caches
result = my_func(3, 4)       # returns cached value
```

Cache files are stored as `{path}/{key}.pickle` where `key` is
`hash_sha256(pickle.dumps((func.__qualname__, args, kwargs)))`.

### `cache_call(*, maxsize=128, get_state=None, is_hash_args=True, path=None, is_sync_node=True, cache=None)` (decorator)

General-purpose caching decorator combining in-memory LRU cache with optional
on-disk persistence.

| Parameter | Default | Description |
|---|---|---|
| `maxsize` | `128` | LRU cache size. Set to `0` to disable in-memory caching |
| `get_state` | `None` | Callable that returns extra state included in the cache key |
| `is_hash_args` | `True` | If `True`, pickle all arguments and use SHA-256 hash as key. If `False`, use `(fname, args, state)` directly (requires `kwargs` to be empty) |
| `path` | `None` | If set, also persist cache to disk at `{path}/data/{key}.pickle` |
| `is_sync_node` | `True` | Disk cache is sync'd across nodes |
| `cache` | `None` | Provide an existing `LRUCache` instance to share across functions |

The decorated function gains two attributes:
- `func.cache` — the underlying `LRUCache` instance.
- `func.clear()` — clears the in-memory cache.

It also accepts `is_force_recompute=True` to bypass the cache.

```python
@q.cache_call(maxsize=128, get_state=q.get_jk_state)
def compute_correlator(cfg_id):
    return heavy_computation(cfg_id)

@q.cache_call(maxsize=64, path="cache/solves", is_sync_node=True)
def solve_propagator(gauge_config, source):
    return invert(gauge_config, source)

# shared cache
shared_cache = q.LRUCache(256)

@q.cache_call(cache=shared_cache)
def func_a(x):
    return x * 2

@q.cache_call(cache=shared_cache)
def func_b(x):
    return x + 1
```

---

## Display Method

### `SetDisplayMethod`

RAII-style class. On construction, sets the display method to
`"py_stdout"` (output via Python's `print`). On destruction (when the object
goes out of scope or is `del`'d), resets to the default display method.

```python
sdm = q.SetDisplayMethod()
q.displayln_info("This prints via Python stdout")
del sdm
```

Equivalent to manual use:

```python
q.set_display_method("py_stdout")
# ... code that prints ...
q.set_display_method()
```

---

## Examples

### Saving and Loading Data

```python
import qlat_utils as q
import numpy as np

data = {
    "config_id": 42,
    "plaquette": np.float32(0.58),
    "correlator": np.array([1.0, 0.8, 0.6, 0.4]),
}

q.save_json_obj(data, "output/config_42.json")
restored = q.load_json_obj("output/config_42.json")
print(restored["plaquette"])    # 0.58 (as numpy.float32)
print(restored["correlator"])   # [1.0, 0.8, 0.6, 0.4] (as numpy.ndarray)
```

### Function Caching on Disk

```python
import qlat_utils as q

@q.pickle_cache("cache/my_results")
def expensive_sum(a, b):
    print("Computing...")
    return a + b

print(expensive_sum(3, 4))   # prints "Computing...", returns 7
print(expensive_sum(3, 4))   # no print, returns cached 7
```

### LRU Caching with State

```python
import qlat_utils as q

counter = 0

def my_state():
    return counter

@q.cache_call(maxsize=128, get_state=my_state)
def compute(x):
    print(f"Computing for x={x} at counter={counter}")
    return x ** 2

print(compute(5))           # computes: "Computing for x=5 at counter=0"
print(compute(5))           # cached, no print
counter = 1
print(compute(5))           # recomputes: "Computing for x=5 at counter=1"
```

### Forced Recompute

```python
@q.cache_call(maxsize=128, path="cache/force_demo")
def expensive(x):
    print("Work...")
    return x * 10

print(expensive(3))                    # computes
print(expensive(3))                    # cached
print(expensive(3, is_force_recompute=True))   # recomputes
```

### Hashing for Cache Keys

```python
import qlat_utils as q
import numpy as np

h1 = q.hash_sha256((32, 64, "gauge"))
h2 = q.hash_sha256(np.array([1.0, 2.0]))
h3 = q.hash_sha256("hello world")
```
