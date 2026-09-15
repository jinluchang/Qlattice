#!/usr/bin/env python3

# Field indexing and NumPy bridge contract tests.
#
# Covers:
#   * Coordinate keys raising TypeError instead of crashing or silently
#     fancy-indexing the site axis,
#   * the flat buffer layout (local site index, first coordinate fastest),
#   * the coordinate-ordered transpose recipe and its inverse,
#   * the writeable zero-copy NumPy view and its lifetime rules.

import numpy as np

import qlat as q

def probe(fn):
    """Return the outcome of calling fn() as a string."""
    try:
        fn()
        return "no-exception"
    except BaseException as e:
        return type(e).__name__


q.begin_with_mpi([q.Coordinate([1, 1, 1, 1])])

total_site = [2, 3, 5, 2]
geo = q.Geometry(q.Coordinate(total_site))
f = q.Field(q.ElemTypeColorMatrix, geo, 2)

### -------------------------------------------------------------------
# Coordinate keys are rejected with a helpful error
### -------------------------------------------------------------------

bad_keys = [
    ("Coordinate", q.Coordinate([0, 0, 0, 0])),
    ("CoordinateD", q.CoordinateD([0.0, 0.0, 0.0, 0.0])),
    ("tuple with Coordinate", (q.Coordinate([1, 1, 0, 0]), 1)),
]
for name, key in bad_keys:
    q.json_results_append(f"getitem {name} = {probe(lambda k=key: f[k])}")

for name, key in bad_keys:
    q.json_results_append(f"setitem {name} = {probe(lambda k=key: f.__setitem__(k, 0))}")

try:
    f[q.Coordinate([0, 0, 0, 0])]
except TypeError as e:
    q.json_results_append(f'message mentions get_elem_xg = {"get_elem_xg" in str(e)}')
    q.json_results_append(
        f'message mentions index_from_coordinate = {"index_from_coordinate" in str(e)}'
    )

# normal NumPy indexing still works, including the slices used in the docs
q.json_results_append(f"f[:] shape = {f[:].shape}")
q.json_results_append(f"f[:] writeable = {f[:].flags.writeable}")
q.json_results_append(f"f[:] equals np.asarray = {np.array_equal(f[:], np.asarray(f))}")
q.json_results_append(f"f[0:2] shape = {f[0:2].shape}")
q.json_results_append(f"f[...,0].shape = {f[..., 0].shape}")
q.json_results_append(f"f[None].shape = {f[None].shape}")
q.json_results_append(f"f[0,0].shape = {f[0, 0].shape}")

### -------------------------------------------------------------------
# buffer layout: flat local site index, first coordinate fastest
### -------------------------------------------------------------------

a = np.asarray(f)
V = a.shape[0]
m = a.shape[1]
elem = a.shape[2:]
q.json_results_append(f"buffer shape = {a.shape}")
q.json_results_append(f"V - local_volume = {V - geo.local_volume}")
q.json_results_append(f"buffer C-contiguous = {a.flags.c_contiguous}")

local_site = geo.node_site.to_list()   # local dimensions (total_site / size_node)
q.json_results_append(f"local site = {local_site}")
q.json_results_append(f"V - prod(local_site) = {V - int(np.prod(local_site))}")

layout_ok = True
for i in range(V):
    xl = geo.coordinate_from_index(i).to_list()
    expect = xl[0] + local_site[0] * (xl[1] + local_site[1] * (xl[2] + local_site[2] * xl[3]))
    if i != expect or geo.index_from_coordinate(q.Coordinate(xl)) != i:
        layout_ok = False
q.json_results_append(f"flat index = coordinate_from_index (first coord fastest) = {layout_ok}")

### -------------------------------------------------------------------
# coordinate-ordered recipe from docs/qlat/qlat_field_indexing.md
### -------------------------------------------------------------------

for i in range(V):
    a[i, 0, 0, 0] = complex(i + 1, 0)

rev = local_site[::-1]
arr = a.reshape(*rev, m, *elem).transpose(3, 2, 1, 0, 4, 5, 6)
q.json_results_append(f"arr axes are local (x,y,z,t,...) = {arr.shape}")

# arr[xl] == buffer at geo.index_from_coordinate(xl), for every local site
coord_ok = all(
    complex(arr[tuple(geo.coordinate_from_index(i).to_list()) + (0, 0, 0)]) == complex(i + 1)
    for i in range(V)
)
q.json_results_append(f"arr[xl] matches flat index = {coord_ok}")

# cross-check against the library's own global-coordinate accessor
xg_list = [geo.coordinate_g_from_l(geo.coordinate_from_index(i)).to_list() for i in range(V)]
gx = np.array(f.get_elem_xg(xg_list, 0)[:]).reshape(V, *elem)
gx_ok = all(
    complex(arr[tuple(geo.coordinate_from_index(i).to_list()) + (0, 0, 0)]) == complex(gx[i, 0, 0])
    for i in range(V)
)
q.json_results_append(f"arr[local coord] == get_elem_xg(global coord) = {gx_ok}")

# inverse recipe round-trips exactly
back = arr.transpose(3, 2, 1, 0, 4, 5, 6).reshape(V, m, *elem)
q.json_results_append(f"inverse recipe round-trips = {np.array_equal(back, a)}")

# write-back through the inverse recipe reaches the field
target = [d // 2 for d in local_site]
cfg = arr.copy()
cfg[tuple(target) + (0, 0, 0)] = complex(999, 0)
a[:] = cfg.transpose(3, 2, 1, 0, 4, 5, 6).reshape(V, m, *elem)
i_target = geo.index_from_coordinate(q.Coordinate(target))
q.json_results_append("write-back lands at index_from_coordinate", complex(a[i_target, 0, 0, 0]) - complex(999, 0))

### -------------------------------------------------------------------
# view semantics and lifetime
### -------------------------------------------------------------------

q.json_results_append(f"view_count before = {f.view_count}")
v1 = f[:]
q.json_results_append(f"view_count with one view = {f.view_count}")
del v1
q.json_results_append(f"view_count after del = {f.view_count}")

v2 = f[:]
q.json_results_append(f"re-init while viewed = {probe(lambda: f.init_from_geo(geo))}")
del v2
q.json_results_append(f"view_count after del 2 = {f.view_count}")

# writing through the view changes the field
w = np.asarray(f)
w[0, 0, 0, 0] = complex(7, 0)
q.json_results_append("write through view reaches field", complex(np.asarray(f)[0, 0, 0, 0]) - complex(7, 0))
del w
q.json_results_append(f"view_count at end = {f.view_count}")

### -------------------------------------------------------------------

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-14)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
