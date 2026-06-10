# `qlat_utils.cache` — Hierarchical In-Memory Cache

Source: `qlat-utils/qlat_utils/cache.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Cache Class](#cache-class)
3. [Global Cache & Management Functions](#global-cache--management-functions)
4. [Examples](#examples)

---

## Overview

The `qlat_utils.cache` module provides a hierarchical, dict-based caching
system built on the `Cache` class (a `dict` subclass). Each `Cache` node
carries a `cache_keys` attribute recording its path from the root, enabling
introspection and targeted cleanup.

A module-level root `cache` object serves as the top of the hierarchy.
Sub-caches are created and retrieved with `mk_cache`. Values can be cleaned
recursively without destroying structure, and the entire Python + C++ cache
stack can be cleared with `clear_all_caches`.

---

## Cache Class

### `class Cache(dict)`

A `dict` subclass that records a hierarchical key path.

| Attribute | Type | Description |
|---|---|---|
| `cache_keys` | `tuple[str, ...]` | Chain of keys from the root cache to this node |

`Cache` behaves exactly like a regular `dict` for item access (`__getitem__`,
`__setitem__`, `pop`, `items`, etc.). The `cache_keys` tuple is purely
informational and used by management functions for logging.

---

## Global Cache & Management Functions

### `cache`

The module-level root `Cache()` instance. All sub-caches are descendants of
this object.

### `mk_cache(*keys, ca=cache)`

Create a nested chain of `Cache` nodes under `ca` for each key in `keys`.
If a sub-cache already exists at any level, it is reused (not recreated).

Returns the innermost `Cache` node.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `*keys` | `str` | — | Hierarchical key path |
| `ca` | `Cache` | `cache` | Parent cache to search/create under |

### `clean_cache(ca=cache)`

Recursively remove all **leaf values** (non-`Cache` entries) from `ca` and its
sub-caches, preserving the `Cache` structure itself. Timed with `@timer`.

### `rm_cache(*keys, ca=cache)`

Remove the sub-cache identified by `keys` from `ca` if it exists. Timed with
`@timer`.

### `list_cache(ca=cache)`

Return a nested `dict` mirroring the cache hierarchy. Leaf values are omitted;
only `Cache` sub-nodes are included. Timed with `@timer`.

### `show_cache_keys(keys)`

Return a human-readable string representation of a `cache_keys` tuple, e.g.
`"['fields_io']['gauge_field']"`.

### `clear_all_caches()`

Clean the Python-level cache (`clean_cache` + `cache.clear()`) and then clear
the C++-level cache via `c.clear_all_caches()`. Timed with `@timer`.

---

## Examples

### Creating and Using a Sub-Cache

```python
import qlat_utils as q

# Create a sub-cache for field I/O handles
cache_fields_io = q.mk_cache("fields_io")

# Store a value keyed by object id
cache_fields_io[id(my_obj)] = (fsel, sbs)

# Retrieve
if id(my_obj) in cache_fields_io:
    c_fsel, c_sbs = cache_fields_io[id(my_obj)]

# Remove
cache_fields_io.pop(id(my_obj), None)
```

### Hierarchical Sub-Caches

```python
import qlat_utils as q

cache_gauge = q.mk_cache("fields_io", "gauge_field")
# Equivalent to:
#   cache["fields_io"]["gauge_field"]
# created step by step, reusing existing nodes.
```

### Cleaning and Removing

```python
import qlat_utils as q

# Remove all values but keep structure
q.clean_cache()

# Remove a specific sub-cache entirely
q.rm_cache("fields_io")  # TODO: rm_cache is buggy -- `for key in keys[-1]` iterates over characters of the string instead of traversing the key path, so it silently returns early without removing anything

# Nuclear option: clear Python + C++ caches
q.clear_all_caches()  # TODO: q.clear_all_caches resolves to the Cython binding from c.py (only clears C++ caches), not the Python version from cache.py that also clears the Python-level cache
```

### Inspecting the Cache

```python
import qlat_utils as q

structure = q.list_cache()
print(structure)  # nested dict of Cache nodes only
```
