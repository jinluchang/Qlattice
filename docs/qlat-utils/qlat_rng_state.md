# `qlat_utils.rng_state` — Deterministic SHA-256-Based Random Number Generator Module

Source: `qlat-utils/qlat_utils/rng_state.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Design Philosophy](#design-philosophy)
3. [RngState Class](#rngstate-class)
   - [Constructor](#constructor)
   - [Copying and Assignment](#copying-and-assignment)
   - [Splitting](#splitting)
   - [Scalar Random Generators](#scalar-random-generators)
   - [Array Random Generators](#array-random-generators)
   - [Selection and Permutation Helpers](#selection-and-permutation-helpers)
4. [Module-Level Utility Functions](#module-level-utility-functions)
5. [C++ Internals Reference](#c-internals-reference)
6. [Examples](#examples)

---

## Overview

`rng_state` is the random number generator module of `qlat_utils`. It provides:

- **`RngState`** — a deterministic, reproducible, split-based PRNG built on
  SHA-256 hashing, designed for parallel lattice-QCD simulations where every
  MPI rank, GPU, or thread must produce independent yet fully reproducible
  random streams without a central coordinator.
- **`get_data_sig()`** — compute a deterministic floating-point signature of
  arbitrary data by dotting with a random ±1 vector.
- **`random_permute()`** — Fisher-Yates shuffle returning a new permuted list.

Key properties of `RngState`:

| Property | Description |
|---|---|
| Determinism | Same seed + same split history = identical sequence on every platform |
| Splitting | Fork child generators by name/index; children are independent of each other |
| No period | SHA-256 based; no short-cycle issues common in LCG/Xorshift generators |
| Thread safe | Each `RngState` instance owns its own state; no global lock needed |
| 64-bit output | `rand_gen()` returns `uint64` values uniformly over `[0, 2^64-1]` |

Python access:

```python
import qlat_utils as q
rs = q.RngState("my-seed")
```

---

## Design Philosophy

### SHA-256 counter mode

Internally, each `RngState` stores a 256-bit SHA-256 hash (`hash[8]` of `uint32`)
plus a byte counter (`numBytes`). When a random number is requested, the
generator increments its `index`, formats a string like `"[index]"` (or
`"[type,index]"` if `type` is set), and feeds it into SHA-256 keyed by the
current hash. The resulting 256-bit digest yields **four** `uint64` values: one
is returned immediately, the other three are cached for subsequent calls. This
gives a 4x throughput boost per SHA-256 invocation.

### Split-based forking

Instead of using a single global seed and skipping ahead (which requires
coordination), `RngState` uses **splitting**: calling `split("label")` produces
a *new* `RngState` whose hash is derived from the parent hash plus the label.
The child's sequence is cryptographically independent of the parent and of all
other children. This is the recommended way to assign independent random streams
to lattice sites, directions, MPI ranks, etc.

### Gaussian via Box-Muller

`g_rand_gen()` uses the Box-Muller transform. One of the two generated Gaussian
values is cached in `gaussian` / `gaussianAvail` and returned on the next call,
so every pair of uniform draws produces two Gaussian draws efficiently.

---

## RngState Class

### Constructor

```python
RngState()
RngState(seed: str | int)
RngState(parent: RngState)
RngState(parent: RngState, seed: str | int)
```

#### Parameters

| Form | Behavior |
|---|---|
| `RngState()` | Create an RNG with the default (zero) initial hash. Useful as a blank slate before assignment. |
| `RngState(seed)` | Create an RNG seeded with `seed`. The seed is converted to a string internally. Equivalent to `RngState()` then `reset(seed)`. |
| `RngState(parent)` | Copy constructor — create an independent but *identical* clone of `parent`. |
| `RngState(parent, seed)` | Split `parent` by `seed` to produce a new independent child. Equivalent to `parent.split(seed)`. |

#### Examples

```python
import qlat_utils as q

rs0 = q.RngState()           # default (zero-hash) state
rs1 = q.RngState("hello")    # seeded with string "hello"
rs2 = q.RngState(42)         # seeded with integer 42 (converted to "42")
rs3 = q.RngState(rs1)        # copy of rs1
rs4 = q.RngState(rs1, "sub") # split child of rs1 with label "sub"
```

---

### Copying and Assignment

#### `copy(is_copying_data=True) -> RngState`

Return a copy of this `RngState`. If `is_copying_data` is `False`, return a
default-initialized (blank) `RngState`.

```python
rs_copy = rs.copy()            # identical copy
rs_blank = rs.copy(False)      # default state, not a copy of rs
```

#### `__copy__` / `__deepcopy__`

Both delegate to `copy()`. Standard `copy.copy(rs)` and `copy.deepcopy(rs)`
work as expected.

#### `__imatmul__` (in-place assignment `@=`)

```python
rs1 @= rs2   # rs1 now has identical state to rs2
```

---

### Splitting

#### `split(seed: str) -> RngState`

Produce a new `RngState` deterministically derived from `self` and `seed`.
The child is fully independent of the parent and of all other children with
different seeds.

```python
rs_root = q.RngState("experiment-1")

rs_site_0   = rs_root.split("site-0")
rs_site_1   = rs_root.split("site-1")
rs_direction = rs_root.split("mu=2")
rs_rank     = rs_root.split("3")
```

**Splitting is the primary mechanism for assigning independent streams in
parallel computations.** Typical pattern for a lattice simulation:

```python
rs_global = q.RngState("simulation-v1")

for site in all_sites:
    rs_site = rs_global.split(site.to_tuple())  # or any unique label
    # use rs_site for all randomness at this site
```

---

### Scalar Random Generators

#### `rand_gen() -> int`

Generate a uniformly distributed random integer in `[0, 2^64 - 1]`.

```python
r = rs.rand_gen()   # e.g. 14873649204851673523
```

#### `u_rand_gen(upper=1.0, lower=0.0) -> float`

Generate a uniformly distributed random `float64` in `[lower, upper)`.

```python
x = rs.u_rand_gen()          # uniform in [0.0, 1.0)
x = rs.u_rand_gen(5.0, 2.0) # uniform in [2.0, 5.0)
```

#### `g_rand_gen(center=0.0, sigma=1.0) -> float`

Generate a Gaussian (normal) distributed random `float64` with the given
`center` (mean) and `sigma` (standard deviation). Uses Box-Muller internally;
caches one of the two generated values for efficiency.

```python
z = rs.g_rand_gen()              # standard normal N(0,1)
x = rs.g_rand_gen(3.0, 0.5)     # N(3.0, 0.5)
```

#### `c_rand_gen(size: Coordinate) -> Coordinate`

Generate a uniformly random lattice coordinate within the hyper-rectangle
defined by `size` (a 4-component `Coordinate`). Useful for picking a random
lattice site given the lattice `total_site`.

```python
total_site = q.Coordinate((16, 8, 8, 8))
coord = rs.c_rand_gen(total_site)  # e.g. Coordinate((7, 3, 5, 2))
```

---

### Array Random Generators

These return NumPy arrays filled with random values of the specified shape.

#### `rand_arr(shape) -> np.ndarray`

Return a `uint64` array of the given `shape` filled with uniform random
integers in `[0, 2^64-1]`.

```python
arr = rs.rand_arr((4, 4))  # shape (4,4), dtype=uint64
```

#### `u_rand_arr(shape) -> np.ndarray`

Return a `float64` array of the given `shape` filled with uniform random
values in `[0.0, 1.0)`.

```python
arr = rs.u_rand_arr((100,))  # 100 uniform floats
```

#### `g_rand_arr(shape) -> np.ndarray`

Return a `float64` array of the given `shape` filled with Gaussian random
values with `center=0.0, sigma=1.0`.

```python
arr = rs.g_rand_arr((8, 8, 8, 8))  # 4D lattice of N(0,1) values
```

---

### Selection and Permutation Helpers

#### `select(l: list) -> Any`

Pick and return a uniformly random element from list `l`.

```python
sites = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
chosen = rs.select(sites)
```

---

## Module-Level Utility Functions

These are free functions in the `qlat_utils` namespace (not methods of
`RngState`).

### `get_data_sig(x, rs: RngState) -> float | complex`

Compute a deterministic signature of data `x` (an ndarray, `LatData`,
`SpinMatrix`, etc.) by dotting the flattened data with a random `{-1, +1}`
vector drawn from `rs`. The result depends only on the *value* of `x`, not its
shape. Useful for checksumming field data.

```python
sig = q.get_data_sig(field_array, rs)
```

### `random_permute(l: list, rs: RngState) -> list`

Return a new list that is a random permutation of `l` (Fisher-Yates shuffle).
Does **not** modify the original list.

```python
permuted = q.random_permute([0, 1, 2, 3, 4], rs)
```

---

## C++ Internals Reference

For developers extending or debugging the C++ layer.

### Struct Layout (`rng-state.h`)

```cpp
struct RngState {
    uint64_t numBytes;       // total bytes processed by SHA-256
    uint32_t hash[8];        // 256-bit SHA-256 state
    uint64_t type;           // optional type tag (ULONG_MAX = unset)
    uint64_t index;          // number of rand_gen calls made
    uint64_t cache[3];       // cached random values (3 unused from last hash)
    RealD    gaussian;       // cached Gaussian value (Box-Muller)
    Int      cacheAvail;     // number of valid entries in cache[]
    bool     gaussianAvail;  // whether `gaussian` holds a valid cached value
};
```

### Free Functions (`rng-state.h`)

| C++ Signature | Description |
|---|---|
| `void reset(RngState& rs)` | Reset to initial (zero-hash) state |
| `void reset(RngState& rs, const std::string& seed)` | Reset and seed from string |
| `void reset(RngState& rs, const Long seed)` | Reset and seed from integer |
| `void split_rng_state(RngState& rs, const RngState& rs0, const std::string& sindex)` | Derive `rs` from `rs0` + label |
| `void split_rng_state(RngState& rs, const RngState& rs0, const Long sindex)` | Same, with integer label |
| `void set_type(RngState& rs, const uint64_t type)` | Set the type tag (must be unset) |
| `uint64_t rand_gen(RngState& rs)` | Return uniform `uint64` |
| `RealD u_rand_gen(RngState& rs, RealD upper, RealD lower)` | Uniform `RealD` in `[lower, upper)` |
| `RealD g_rand_gen(RngState& rs, RealD center, RealD sigma)` | Gaussian `RealD` |
| `void compute_hash_with_input(uint32_t hash[8], const RngState& rs, const std::string& input)` | Core SHA-256 computation |

### Member Functions

| Method | Description |
|---|---|
| `rs.split(sindex)` | Return new `RngState` split by string label |
| `rs.newtype(type)` | Return a copy with the type tag set |

### How `rand_gen` Works Internally

1. Increment `index`.
2. If `cacheAvail > 0`, return `cache[--cacheAvail]`.
3. Otherwise, compute `SHA-256(hash ‖ "[index]")` (or `"[type,index]"` if `type` is set).
4. The 256-bit digest yields 4 × `uint64`. Return the last one immediately; cache the first three.
5. Update `hash` and increment `numBytes`.

### How `split_rng_state` Works Internally

1. Format a label string: `"[index] {sindex}"` or `"[type,index] {sindex}"`.
2. Pad the label to a multiple of 64 bytes.
3. Feed the padded label through `SHA-256::processBlock` keyed by `rs0.hash`.
4. The new `rs` gets the resulting hash, `numBytes += nBlocks * 64`, and all
   cache/gaussian state is cleared.

---

## Examples

### Basic Usage

```python
import qlat_utils as q

rs = q.RngState("my-seed")

# Single random values
print(rs.rand_gen())      # random uint64
print(rs.u_rand_gen())    # random float in [0, 1)
print(rs.g_rand_gen())    # random Gaussian N(0,1)

# Arrays
arr = rs.u_rand_arr((3, 3))  # 3x3 uniform float array
print(arr)
```

### Splitting for Parallel Sites

```python
import qlat_utils as q

rs = q.RngState("lattice-run-1")

# Assign independent streams per site
total_site = q.Coordinate((8, 8, 8, 8))
for t in range(total_site[0]):
    for z in range(total_site[1]):
        for y in range(total_site[2]):
            for x in range(total_site[3]):
                rs_site = rs.split(f"({t},{z},{y},{x})")
                # Use rs_site for all randomness at this site
                link_value = rs_site.g_rand_gen()
```

### Reproducible Sub-streams

```python
import qlat_utils as q

rs = q.RngState("experiment")

# These two calls always produce the same results regardless of ordering
rs_a = rs.split("stream-a")
rs_b = rs.split("stream-b")

# Each stream is deterministic
val_a = rs_a.u_rand_gen()  # always the same
val_b = rs_b.u_rand_gen()  # always the same, independent of val_a
```

### Random Permutation

```python
import qlat_utils as q

rs = q.RngState("shuffle-seed")
sites = list(range(64))
permuted_sites = q.random_permute(sites, rs)
print(permuted_sites)  # deterministic permutation
```

### Checksumming with `get_data_sig`

```python
import qlat_utils as q
import numpy as np

rs = q.RngState("checksum")
data = np.random.randn(100)
sig = q.get_data_sig(data, rs)
# sig is a deterministic float depending on data values
```
