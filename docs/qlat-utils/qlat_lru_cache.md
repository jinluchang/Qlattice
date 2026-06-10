# `qlat_utils.lru_cache` — Fixed-Capacity LRU Cache

Source: `qlat-utils/qlat_utils/lru_cache.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [API](#api)
3. [Examples](#examples)

---

## Overview

The `qlat_utils.lru_cache` module provides `LRUCache`, a fixed-capacity
Least Recently Used cache backed by `collections.OrderedDict`. All operations
(lookup, insertion, update) run in O(1) time.

When the cache reaches its capacity, inserting a new item automatically evicts
the least recently accessed entry. Accessing an item (via `get` or `[]`) moves
it to the most-recently-used position.

---

## API

### `class LRUCache(capacity: int)`

| Parameter | Type | Description |
|---|---|---|
| `capacity` | `int` | Maximum number of entries the cache can hold |

#### Attributes

| Attribute | Type | Description |
|---|---|---|
| `cache` | `OrderedDict` | Internal storage |
| `capacity` | `int` | Maximum number of entries |

#### `get(key, default=None)`

Return the value for `key`, or `default` if not found. Moves the key to the
most-recently-used position on hit.

#### `__getitem__(key)`

Return the value for `key`. Raises `KeyError` if not found. Moves the key to
the most-recently-used position.

#### `__setitem__(key, value)`

Insert or update `key` with `value`. Moves the key to the most-recently-used
position. If the cache exceeds `capacity`, evicts the least recently used
entry.

#### `__contains__(key)`

Return `True` if `key` is in the cache. Does **not** affect recency order.

#### `clear()`

Remove all entries from the cache.

---

## Examples

### Basic Usage

```python
import qlat_utils as q

cache = q.LRUCache(capacity=3)

cache["a"] = 1
cache["b"] = 2
cache["c"] = 3

print(cache["a"])  # 1; "a" is now most recent

cache["d"] = 4     # evicts "b" (least recent)

print("b" in cache)  # False
print(cache.get("b"))  # None
```

### Gauge Field Caching

```python
import qlat_utils as q

# Cache up to 4 gauge configurations
gauge_cache = q.LRUCache(capacity=4)

for traj in range(10):
    gf = load_gauge_field(traj)  # hypothetical loader
    gauge_cache[traj] = gf
    # Old trajectories are evicted automatically
```

### Checking Eviction Order

```python
import qlat_utils as q

c = q.LRUCache(2)
c["x"] = 10
c["y"] = 20
_ = c["x"]       # touch "x" → "y" is now LRU
c["z"] = 30      # evicts "y"

print("y" in c)  # False
print(c["x"])    # 10
print(c["z"])    # 30
```
