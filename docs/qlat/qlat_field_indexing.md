# `qlat` field indexing and the NumPy bridge

Source: `qlat/qlat/field_base.pyx`, `qlat/qlat/field_types.pyx`

> **Note:** Update this document when updating the source files.

This page covers the three indexing conventions that trip people up when moving
between `qlat` fields and NumPy arrays:

1. the **flat local site index** used by `f[i]`,
2. the **global coordinate** (`_xg`) API,
3. the **coordinate-ordered** array layout for analysis code.

It also covers initialisation traps and the view/assign idioms.

---

## 1. Flat local site index

`np.asarray(f)` (equivalently `f[:]`) returns a writeable, zero-copy view with
shape

```
(local_volume, multiplicity, *elem_shape)
```

C-contiguous. Axis 0 is the **flat local site index**: the first coordinate
varies fastest, matching `geo.coordinate_from_index`. The site axis covers the
**local** volume on each MPI rank, i.e. `geo.node_site`, not the global
`total_site`. For local dimensions `node_site = [Lx, Ly, Lz, Lt]`, local
coordinate `(x, y, z, t)` maps to

```
i = x + Lx * (y + Ly * (z + Lz * t))
```

With a single rank, `node_site == total_site`; with several ranks each rank's
buffer has shape `(prod(node_site), multiplicity, *elem_shape)`.

```python
import numpy as np
import qlat as q

q.begin_with_mpi()
geo = q.Geometry(q.Coordinate([4, 4, 4, 4]))
f = q.Field(q.ElemTypeColorMatrix, geo, 1)

i = geo.index_from_coordinate(q.Coordinate([1, 2, 3, 0]))
xl = geo.coordinate_from_index(i)          # Coordinate([1, 2, 3, 0])

a = np.asarray(f)                          # (local_volume, 1, 3, 3)
a[i, 0]                                    # the element at local (1, 2, 3, 0)
f[i, 0] = q.ColorMatrix()                  # write through
```

`f[i]` and `f[i, m]` return **views**, so `f[i][0, 0] = 1.0` writes to the
field. A `Coordinate` is **not** an accepted index — passing one raises
`TypeError` with a pointer to the APIs below.

---

## 2. Global coordinates (`_xg`)

The `_xg` methods take **global** coordinates (identical on every MPI rank) and
are collective:

| method | returns |
|---|---|
| `f.get_elem_xg(xg_arr, m)` | `SelectedPoints` with shape `(N, 1, *elem_shape)` |
| `f.set_elem_xg(xg_arr, m, val)` | writes element `m` at the coordinates |
| `f.get_elems_xg(xg_arr)` | `SelectedPoints` with all multiplicities |
| `f.set_elems_xg(xg_arr, val)` | writes all multiplicities |

```python
sp = f.get_elem_xg([[0, 0, 0, 0], [1, 2, 3, 0]], 0)
sp.shape        # (2, 1, 3, 3)
val = sp[:]     # (2, 1, 3, 3) NumPy array
```

> **Collective:** `xg_arr` must be identical on all MPI processes. Passing
> different values on different ranks produces undefined results.

On a single rank with no decomposition, a global coordinate equals a local one;
with multiple ranks, `geo.coordinate_l_from_g(xg)` converts.

---

## 3. Coordinate-ordered arrays

The buffer is flat over sites, so `a.reshape(Lx, Ly, Lz, Lt, m, ...)` does
**not** give an `[x, y, z, t]`-indexed array: the flat site axis already has the
first coordinate varying fastest, which is the opposite of NumPy's
C-contiguous convention for a `(Lx, Ly, Lz, Lt, ...)` array. The reshape must
therefore use the reversed lattice dimensions, and then the four space-time
axes are reversed back. Use the **local** dimensions from `geo.node_site`, not
the global `total_site`:

```python
loc = geo.node_site.to_list()           # [Lx, Ly, Lz, Lt] on this rank
rev = loc[::-1]                         # [Lt, Lz, Ly, Lx]
a = np.asarray(f)                       # (V, m, *elem_shape), V = prod(loc)
elem = a.shape[2:]                      # e.g. (3, 3)

# flat buffer -> arr[x, y, z, t, m, *elem]
arr = a.reshape(*rev, a.shape[1], *elem).transpose(3, 2, 1, 0, 4, 5, 6)

# sanity check: arr index is the geometric local coordinate
xl = geo.coordinate_from_index(0)
assert arr[tuple(xl.to_list()) + (0, 0, 0)] == f[geo.index_from_coordinate(xl), 0, 0, 0]
```

The result is a view with axes `(Lx, Ly, Lz, Lt, multiplicity, *elem_shape)` —
so `arr[x, y, z, t]` is the element at that local coordinate. To cross-check
against the library, compare with
`f.get_elem_xg([geo.coordinate_g_from_l(xl).to_list()], m)`, which takes
**global** coordinates.

### Writing back

The inverse operation needs **no transpose after** the reshape:

```python
cfg = arr.copy()                        # edit at will, indexed [x, y, z, t]
idx = [d // 2 for d in loc]
cfg[tuple(idx) + (0, 0, 0)] = q.ColorMatrix()

a[:] = cfg.transpose(3, 2, 1, 0, 4, 5, 6).reshape(a.shape)
```

The transpose appears only on the coordinate-ordered side; it is its own
inverse for the four space-time axes.

For a pure copy (rather than a view), use `.copy()` on `arr`.

---

## 4. Initialisation

`qlat` needs its global geometry node before a `Geometry` can be built from a
`total_site`, so:

```python
import qlat as q

q.begin_with_mpi()                       # or q.begin(id_node, size_node)
geo = q.Geometry(q.Coordinate([4, 4, 4, 4]))
```

Common traps:

- `q.Geometry(total_site)` before `begin*` raises `RuntimeError` (previously it
  crashed with `SIGFPE`).
- `mpirun -n 1 python3 script.py` does **not** help: `mpirun` only sets the
  environment, the process still has to call `MPI_Init`, and `qlat`
  additionally needs `begin_comm` to populate `geon`.
- `Geometry(id_node, size_node, node_site)` works without `begin*` because it
  does not touch the global `geon` — useful for testing.
- Calling `begin_with_mpi()` twice pushes a second communicator; call it once.

---

## 5. Summary of conventions

| API | index meaning | layout |
|---|---|---|
| `f[i]`, `f[i, m]`, `np.asarray(f)` | flat local site index | `(local_volume, multiplicity, *elem_shape)`, first coordinate fastest |
| `f.get_elem_xg(xg, m)` | global coordinate | returns `SelectedPoints`, shape `(N, 1, *elem_shape)` |
| coordinate-ordered `arr` (section 3) | `[x, y, z, t]` | `(Lx, Ly, Lz, Lt, multiplicity, *elem_shape)` |
| `SelectedField` / `SelectedPoints` buffer | selection index | `(n_elems, multiplicity, *elem_shape)` — axis 0 is a selection index, not a lattice coordinate |
