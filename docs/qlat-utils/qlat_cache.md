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
my_obj = object()
cache_fields_io[id(my_obj)] = ("fsel_data", "sbs_data")

# Retrieve
assert id(my_obj) in cache_fields_io
c_fsel, c_sbs = cache_fields_io[id(my_obj)]
assert c_fsel == "fsel_data"
assert c_sbs == "sbs_data"

# Remove
cache_fields_io.pop(id(my_obj), None)
assert id(my_obj) not in cache_fields_io
```

### Hierarchical Sub-Caches

```python
import qlat_utils as q

cache_gauge = q.mk_cache("fields_io", "gauge_field")
# Equivalent to:
#   cache["fields_io"]["gauge_field"]
# created step by step, reusing existing nodes.

# Store and retrieve a value
cache_gauge["config_100"] = {"beta": 6.0, "L": 8}
assert cache_gauge["config_100"]["beta"] == 6.0

# mk_cache reuses existing nodes (does not overwrite)
cache_gauge2 = q.mk_cache("fields_io", "gauge_field")
assert cache_gauge is cache_gauge2
```

### Removing a Sub-Cache

```python
import qlat_utils as q

# Set up a nested cache
q.mk_cache("fields_io", "gauge_field")
q.mk_cache("fields_io", "prop_field")

# Remove one sub-cache (keeps sibling)
q.rm_cache("fields_io", "gauge_field")

# Verify it is gone
assert "gauge_field" not in q.cache.get("fields_io", {})
# Sibling still exists
assert "prop_field" in q.cache["fields_io"]
```

### Cleaning Values (Preserving Structure)

```python
import qlat_utils as q

# Set up cache with values
cache_io = q.mk_cache("fields_io")
cache_io["key_a"] = "value_a"
cache_io["key_b"] = "value_b"

# clean_cache removes leaf values but keeps Cache sub-nodes
q.clean_cache()
assert len(cache_io) == 0
# The Cache node itself still exists in the hierarchy
assert "fields_io" in q.cache
```

### Clearing All Caches

```python
import qlat_utils as q

# Set up some caches with values
cache_io = q.mk_cache("fields_io")
cache_io["key"] = "value"

# clear_all_caches clears both Python-level and C++-level caches
q.clear_all_caches()

# Python cache is empty after clearing
assert len(q.cache) == 0
```

### Inspecting the Cache

```python
import qlat_utils as q

# Create some nested caches
q.mk_cache("fields_io", "gauge_field")
q.mk_cache("fields_io", "prop_field")

structure = q.list_cache()
# Returns a nested dict with Cache sub-nodes only (no leaf values)
assert "fields_io" in structure
assert "gauge_field" in structure["fields_io"]
assert "prop_field" in structure["fields_io"]
```

### Multiple Independent Cache Trees

```python
import qlat_utils as q

# You can create independent cache trees under different roots
root_a = q.Cache()
root_b = q.Cache()

ca = q.mk_cache("group", "sub", ca=root_a)
cb = q.mk_cache("group", "sub", ca=root_b)

ca["from_a"] = 1
cb["from_b"] = 2

assert ca["from_a"] == 1
assert cb["from_b"] == 2
# They are independent
assert "from_b" not in ca
assert "from_a" not in cb
```

### Cache Keys Introspection

```python
import qlat_utils as q

cache_gauge = q.mk_cache("fields_io", "gauge_field")
# cache_keys records the path from the root
assert cache_gauge.cache_keys == ("fields_io", "gauge_field")

# show_cache_keys produces a human-readable string
print(q.show_cache_keys(cache_gauge.cache_keys))
# Output: ['fields_io']['gauge_field']
```
