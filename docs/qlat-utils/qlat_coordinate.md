# `qlat_utils.coordinate` — 4-Component Coordinate Types for Lattice QCD

Source: `qlat-utils/qlat_utils/coordinate.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Coordinate Class](#coordinate-class)
   - [Constructors](#constructors)
   - [Copying and Assignment](#copying-and-assignment)
   - [Serialization](#serialization)
   - [Norms and Volumes](#norms-and-volumes)
   - [Index Conversion](#index-conversion)
   - [Arithmetic](#arithmetic)
   - [Indexing and Iteration](#indexing-and-iteration)
   - [Comparison](#comparison)
3. [CoordinateD Class](#coordinated-class)
   - [Constructors](#constructors-1)
   - [Serialization](#serialization-1)
   - [Norm](#norm)
   - [Arithmetic](#arithmetic-1)
   - [Indexing and Comparison](#indexing-and-comparison)
4. [Module-Level Functions](#module-level-functions)
   - [Integer Functions](#integer-functions)
   - [Double Functions](#double-functions)
5. [Examples](#examples)

---

## Overview

The `qlat_utils.coordinate` module provides two 4-component coordinate types for
representing positions and extents on a 4-dimensional lattice:

- **`Coordinate`** — integer coordinates (`int32`), used for lattice site
  positions, lattice sizes, and block dimensions.
- **`CoordinateD`** — double-precision coordinates (`float64`), used for
  continuous momenta, fractional positions, and interpolation points.

Both types support element-wise arithmetic, modular operations, serialization,
and iteration. The integer `Coordinate` additionally provides index↔coordinate
conversion for flat-array access and spatial volume computation.

```python
import qlat_utils as q

total_site = q.Coordinate([8, 8, 8, 8])
site = q.Coordinate([3, 5, 2, 7])
mom = q.CoordinateD([1.0, 2.0, 3.0, 0.5])
```

---

## Coordinate Class

`Coordinate` is a 4-component integer coordinate. Internally it wraps a
C++ `array<Int, 4>` where `Int` is `int32_t`. Components are indexed 0–3,
conventionally representing `(x, y, z, t)` or `(μ=0,1,2,3)`.

### Constructors

```python
Coordinate()
Coordinate(x: list | tuple | np.ndarray | Coordinate)
```

| Form | Behavior |
|---|---|
| `Coordinate()` | Zero-initialize all components. |
| `Coordinate([3, 5, 2, 7])` | From a list of exactly 4 integers. |
| `Coordinate((3, 5, 2, 7))` | From a tuple of exactly 4 integers. |
| `Coordinate(np.array([3, 5, 2, 7]))` | From a shape-(4,) numpy array. |
| `Coordinate(other)` | Copy constructor from another `Coordinate`. |

```python
import qlat_utils as q

c0 = q.Coordinate()                # (0, 0, 0, 0)
c1 = q.Coordinate([3, 5, 2, 7])    # from list
c2 = q.Coordinate((8, 8, 8, 8))    # from tuple
c3 = q.Coordinate(c1)               # copy of c1
```

### Copying and Assignment

```python
c.copy(is_copying_data=True) -> Coordinate
c @= other                         # in-place assignment
```

| Method | Description |
|---|---|
| `copy(is_copying_data=True)` | Return a copy. If `is_copying_data` is `False`, return a zero-initialized `Coordinate`. |
| `__copy__()` | Delegates to `copy()`. |
| `__deepcopy__(memo)` | Delegates to `copy()`. |
| `__imatmul__(other)` | `c @= other` — copy the value of `other` into `c`. |

```python
c1 = q.Coordinate([1, 2, 3, 4])
c2 = c1.copy()                     # independent copy
c2[0] = 99
assert c1[0] == 1                  # unaffected

c1 @= q.Coordinate([5, 6, 7, 8])   # in-place assignment
```

### Serialization

| Method | Returns | Description |
|---|---|---|
| `to_list()` | `list[int]` | Components as a Python list. |
| `to_tuple()` | `tuple[int]` | Components as a Python tuple. |
| `to_numpy()` | `np.ndarray` | Components as a shape-(4,) `int32` array. |
| `from_list(x)` | `None` | Set components from list/tuple/ndarray of length 4. |

`Coordinate` also supports `__getstate__` / `__setstate__` for Python pickling,
and `__repr__` for readable display.

```python
c = q.Coordinate([3, 5, 2, 7])

assert c.to_list() == [3, 5, 2, 7]
assert c.to_tuple() == (3, 5, 2, 7)
assert np.array_equal(c.to_numpy(), np.array([3, 5, 2, 7], dtype=np.int32))

c.from_list([99, 88, 77, 66])
assert c.to_list() == [99, 88, 77, 66]
```

### Norms and Volumes

| Method | Returns | Description |
|---|---|---|
| `sqr()` | `Long` (`int64`) | Sum of squares: `x[0]² + x[1]² + x[2]² + x[3]²`. |
| `r_sqr()` | `int` | Spatial distance squared: `x[0]² + x[1]² + x[2]²`. |
| `volume()` | `Long` (`int64`) | Product of all four components. |
| `spatial_volume()` | `Long` (`int64`) | Product of first three components (spatial only). |

```python
lattice = q.Coordinate([8, 8, 8, 16])

assert lattice.sqr() == 8*8 + 8*8 + 8*8 + 16*16
assert lattice.r_sqr() == 8*8 + 8*8 + 8*8
assert lattice.volume() == 8 * 8 * 8 * 16
assert lattice.spatial_volume() == 8 * 8 * 8
```

### Index Conversion

Map between 4D coordinates and 1D flat indices (row-major with component 0
fastest-varying):

```python
c.from_index(index: int, size: Coordinate)
c.to_index(size: Coordinate) -> int
```

| Method | Description |
|---|---|
| `from_index(index, size)` | Set `c` to the coordinate corresponding to flat `index` within bounds `[0, size)`. |
| `to_index(size)` | Return the flat index of `c` given box size `size`. |

```python
size = q.Coordinate([4, 4, 4, 8])

c = q.Coordinate()
c.from_index(0, size)
assert c.to_list() == [0, 0, 0, 0]

c.from_index(1, size)
assert c.to_list() == [1, 0, 0, 0]

c.from_index(4, size)
assert c.to_list() == [0, 1, 0, 0]

# Round-trip
c = q.Coordinate([3, 2, 1, 5])
idx = c.to_index(size)
c2 = q.Coordinate()
c2.from_index(idx, size)
assert c == c2
```

### Arithmetic

All arithmetic operators work element-wise on the four components.

| Operator | Description |
|---|---|
| `c1 + c2` | Element-wise addition. |
| `c1 - c2` | Element-wise subtraction. |
| `-c` | Unary negation. |
| `+c` | Unary plus (identity). |
| `c * n` / `n * c` | Scalar multiplication (`n` is an `int`). |
| `c1 * c2` | Element-wise (Hadamard) product. |
| `c1 % c2` | Element-wise modulo. |
| `c1 // c2` | Element-wise integer division (floor). |
| `c1 / n` (C++) | Element-wise scalar division (not exposed in Python). |

```python
a = q.Coordinate([1, 2, 3, 4])
b = q.Coordinate([5, 6, 7, 8])

assert (a + b).to_list() == [6, 8, 10, 12]
assert (b - a).to_list() == [4, 4, 4, 4]
assert (-a).to_list() == [-1, -2, -3, -4]
assert (a * 3).to_list() == [3, 6, 9, 12]
assert (a * b).to_list() == [5, 12, 21, 32]
assert (b % a).to_list() == [0, 0, 1, 0]
assert (b // q.Coordinate([2, 2, 2, 2])).to_list() == [2, 3, 3, 4]
```

### Indexing and Iteration

```python
c[k]           # get component k (0 ≤ k < 4)
c[k] = val     # set component k
for x in c:    # iterate over 4 components
```

```python
c = q.Coordinate([3, 5, 2, 7])
assert c[0] == 3
assert c[3] == 7

c[2] = 99
assert c[2] == 99

vals = list(c)
assert vals == [3, 5, 99, 7]
```

### Comparison

```python
c1 == c2     # element-wise equality
```

Only `==` is supported. `!=`, `<`, `>`, etc. are not exposed in Python.

```python
a = q.Coordinate([1, 2, 3, 4])
b = q.Coordinate([1, 2, 3, 4])
c = q.Coordinate([0, 0, 0, 0])

assert a == b
assert not (a == c)
```

---

## CoordinateD Class

`CoordinateD` is the double-precision counterpart of `Coordinate`. It wraps a
C++ `array<RealD, 4>` where `RealD` is `double`. It shares the same serialization
API as `Coordinate` but uses `float64` instead of `int32` and has no
index-conversion or volume-computation methods.

### Constructors

```python
CoordinateD()
CoordinateD(x: list | tuple | np.ndarray | CoordinateD | Coordinate)
```

| Form | Behavior |
|---|---|
| `CoordinateD()` | Zero-initialize all components. |
| `CoordinateD([1.0, 2.0, 3.0, 0.5])` | From a list of 4 numbers. |
| `CoordinateD(c)` | Copy constructor from `CoordinateD` or `Coordinate` (integers are promoted to doubles). |

```python
import qlat_utils as q

d0 = q.CoordinateD()                       # (0.0, 0.0, 0.0, 0.0)
d1 = q.CoordinateD([1.0, 2.0, 3.0, 0.5])   # from list
d2 = q.CoordinateD((4.0, 5.0, 6.0, 7.0))    # from tuple
d3 = q.CoordinateD(d1)                       # copy of d1

# Convert integer Coordinate to CoordinateD
c = q.Coordinate([3, 5, 2, 7])
d4 = q.CoordinateD(c)                        # (3.0, 5.0, 2.0, 7.0)
```

### Serialization

Same API as `Coordinate` but returns `float64` values:

```python
d = q.CoordinateD([1.0, 2.0, 3.0, 0.5])

assert d.to_list() == [1.0, 2.0, 3.0, 0.5]
assert d.to_tuple() == (1.0, 2.0, 3.0, 0.5)
assert np.array_equal(d.to_numpy(), np.array([1.0, 2.0, 3.0, 0.5], dtype=np.float64))

d.from_list([99.0, 88.0, 77.0, 66.0])
assert d.to_list() == [99.0, 88.0, 77.0, 66.0]
```

### Norm

| Method | Returns | Description |
|---|---|---|
| `sqr()` | `RealD` | Sum of squares: `x[0]² + x[1]² + x[2]² + x[3]²`. |

`CoordinateD` does **not** have `r_sqr()`, `volume()`, or `spatial_volume()`.

```python
d = q.CoordinateD([3.0, 4.0, 0.0, 0.0])
assert abs(d.sqr() - 25.0) < 1e-15
```

### Arithmetic

| Operator | Description |
|---|---|
| `d1 + d2` | Element-wise addition. |
| `d1 - d2` | Element-wise subtraction. |
| `-d` | Unary negation. |
| `d * n` / `n * d` | Scalar multiplication (`n` is `int` or `float`). |
| `d1 * d2` | Element-wise (Hadamard) product. |
| `d1 / d2` | True element-wise division (`__truediv__`). |
| `d1 % d2` | Element-wise modulo. |

```python
a = q.CoordinateD([1.0, 2.0, 3.0, 4.0])
b = q.CoordinateD([5.0, 6.0, 7.0, 8.0])

assert (a + b).to_list() == [6.0, 8.0, 10.0, 12.0]
assert (a * 0.5).to_list() == [0.5, 1.0, 1.5, 2.0]
assert (a / q.CoordinateD([2.0, 2.0, 2.0, 2.0])).to_list() == [0.5, 1.0, 1.5, 2.0]
```

### Indexing and Comparison

Same as `Coordinate`: `d[k]`, `d[k] = val`, `for x in d`, `d1 == d2`.

```python
d = q.CoordinateD([1.0, 2.0, 3.0, 4.0])
d[0] = 9.0
assert d[0] == 9.0
assert d == q.CoordinateD([9.0, 2.0, 3.0, 4.0])
```

---

## Module-Level Functions

### Integer Functions

#### `mod_coordinate(c, size) -> Coordinate`

Wrap each component into `[0, size[i])` via element-wise modulo. Equivalent to
`c % size`.

```python
c = q.Coordinate([9, -3, 15, 5])
size = q.Coordinate([4, 4, 4, 8])
result = q.mod_coordinate(c, size)
assert result.to_list() == [1, 1, 3, 5]
```

#### `smod_coordinate(c, size) -> Coordinate`

Shortest-distance modulo: wrap each component into `[-size[i]/2, size[i]/2)`.

```python
c = q.Coordinate([3, 0, -2, 7])
size = q.Coordinate([4, 4, 4, 8])
result = q.smod_coordinate(c, size)
assert result.to_list() == [-1, 0, -2, -1]
```

#### `smod_sym_coordinate(c, size) -> Coordinate`

Symmetric `smod`: like `smod` but components equal to `size[i]/2` become `0`.

```python
c = q.Coordinate([2, 2, 0, 4])
size = q.Coordinate([4, 4, 4, 8])
result = q.smod_sym_coordinate(c, size)
assert result.to_list() == [0, 0, 0, 0]
```

#### `middle_mod_coordinate(x, y, size) -> Coordinate`

Return the midpoint of two coordinates on a periodic lattice: wraps both into
`[0, size)`, computes the shortest-distance midpoint, and wraps the result back
into `[0, size)`.

```python
size = q.Coordinate([8, 8, 8, 8])
a = q.Coordinate([1, 1, 1, 1])
b = q.Coordinate([5, 5, 5, 5])
mid = q.middle_mod_coordinate(a, b, size)
# Midpoint accounting for periodicity
```

#### `coordinate_from_index(index, size) -> Coordinate`

Convert a flat index to a `Coordinate` within bounds `[0, size)`.

```python
size = q.Coordinate([4, 4, 4, 8])
c = q.coordinate_from_index(7, size)
assert c.to_list() == [3, 1, 0, 0]
```

#### `index_from_coordinate(x, size) -> int`

Convert a `Coordinate` to a flat index. The coordinate can be negative; the
function applies `mod(x, size)` internally.

```python
size = q.Coordinate([4, 4, 4, 8])
c = q.Coordinate([3, 1, 0, 0])
idx = q.index_from_coordinate(c, size)
assert idx == 7
assert q.coordinate_from_index(idx, size) == c
```

### Double Functions

#### `mod_coordinate_d(c, size) -> CoordinateD`

Double-precision element-wise modulo. Equivalent to `c % size`.

```python
c = q.CoordinateD([9.5, -1.0, 0.0, 5.0])
size = q.CoordinateD([4.0, 4.0, 4.0, 8.0])
result = q.mod_coordinate_d(c, size)
```

#### `smod_coordinate_d(c, size) -> CoordinateD`

Shortest-distance modulo for double coordinates.

```python
c = q.CoordinateD([3.5, 0.0, -2.0, 7.0])
size = q.CoordinateD([4.0, 4.0, 4.0, 8.0])
result = q.smod_coordinate_d(c, size)
```

#### `smod_sym_coordinate_d(c, size) -> CoordinateD`

Symmetric `smod` for double coordinates.

```python
c = q.CoordinateD([2.0, 2.0, 0.0, 4.0])
size = q.CoordinateD([4.0, 4.0, 4.0, 8.0])
result = q.smod_sym_coordinate_d(c, size)
```

#### `middle_mod_coordinate_d(x, y, size) -> CoordinateD`

Double-precision midpoint calculation on a periodic lattice.

```python
size = q.CoordinateD([8.0, 8.0, 8.0, 8.0])
a = q.CoordinateD([1.0, 1.0, 1.0, 1.0])
b = q.CoordinateD([5.0, 5.0, 5.0, 5.0])
mid = q.middle_mod_coordinate_d(a, b, size)
```

---

## Examples

### Basic Usage

```python
import qlat_utils as q

size = q.Coordinate([8, 8, 8, 16])
site = q.Coordinate([3, 5, 2, 7])

print(f"Lattice size: {size.to_list()}")
print(f"Site: {site}")

print(f"Volume: {size.volume()}")
print(f"Spatial volume: {size.spatial_volume()}")
```

### Iterating Over a Lattice

```python
import qlat_utils as q

size = q.Coordinate([4, 4, 4, 8])
for idx in range(size.volume()):
    c = q.coordinate_from_index(idx, size)
    assert q.index_from_coordinate(c, size) == idx
```

### Coordinate Arithmetic

```python
import qlat_utils as q

size = q.Coordinate([8, 8, 8, 8])
site = q.Coordinate([3, 5, 2, 7])
hop = q.Coordinate([1, 0, 0, 0])   # hop in +x direction

neighbor = site + hop
assert neighbor == q.Coordinate([4, 5, 2, 7])

opposite = -hop
assert opposite == q.Coordinate([-1, 0, 0, 0])
```

### Modular Arithmetic on a Periodic Lattice

```python
import qlat_utils as q

size = q.Coordinate([8, 8, 8, 8])

# Wrap a coordinate into fundamental domain
c = q.Coordinate([9, -3, 15, 7])
wrapped = q.mod_coordinate(c, size)
assert wrapped == q.Coordinate([1, 5, 7, 7])

# Shortest-distance displacement
disp = q.smod_coordinate(q.Coordinate([7, 0, 2, 3]), size)
assert disp == q.Coordinate([-1, 0, 2, 3])
```

### Round-Trip Index Conversion

```python
import qlat_utils as q

size = q.Coordinate([4, 4, 4, 8])

for i in [0, 1, 2, 3, 4, 7, 8, 15, 63, 64, 255, 511]:
    c = q.coordinate_from_index(i, size)
    assert q.index_from_coordinate(c, size) == i
```

### Using CoordinateD for Momenta

```python
import qlat_utils as q

size = q.Coordinate([8, 8, 8, 16])
mom = q.CoordinateD([1.0, 2.0, 3.0, 0.0])
mom_size = q.CoordinateD([8.0, 8.0, 8.0, 16.0])

# Scale momentum by 2π/L
phase = mom * (2.0 * 3.141592653589793)
phase = phase / mom_size

print(f"Momentum: {mom.to_list()}")
print(f"Phase: {phase.to_list()}")
```
