#!/usr/bin/env python3

# Coordinate / CoordinateD sequence protocol tests.
#
# Covers the fixes for:
#   * c[-1]    previously OverflowError  -> now c[3]
#   * c[4]     previously AssertionError -> now IndexError (and not stripped by -O)
#   * c[1.5]   previously silently c[1]  -> now TypeError
#   * c[1] = 1.5 previously silently truncated to 1 -> now TypeError
#   * c[1:3]   previously TypeError      -> still TypeError (integers only)

import numpy as np

import qlat as q


def jr_int(label, *vals):
    """json_results_append with integer values (which must be passed as floats)."""
    q.json_results_append(label, np.array([float(v) for v in vals]))


# probes must record numeric values; encode the outcome instead of a string
OUTCOME_CODE = {
    "IndexError": 1.0,
    "TypeError": 2.0,
    "OverflowError": 3.0,
    "AssertionError": 4.0,
    "no-exception": 0.0,
}


def probe(fn):
    """Return a numeric code for the outcome of calling fn()."""
    try:
        fn()
        return OUTCOME_CODE["no-exception"]
    except BaseException as e:
        return OUTCOME_CODE.get(type(e).__name__, 9.0)

q.begin_with_mpi()

c = q.Coordinate([1, 2, 3, 4])
d = q.CoordinateD([1.5, 2.5, 3.5, 4.5])

jr_int("Coordinate component access", c[0], c[1], c[2], c[3])

# negative indexing
jr_int("Coordinate negative index", c[-1], c[-2], c[-3], c[-4])

# numpy integer keys (__index__)
jr_int("Coordinate numpy int key", int(c[np.int64(2)]))

# out of range
for key in (4, 5, -5, 100):
    q.json_results_append(f"Coordinate c[{key}]", probe(lambda k=key: c[k]))

# non-integers
for key in (1.5, 1.0, "1", None):
    q.json_results_append(f"Coordinate c[{key!r}]", probe(lambda k=key: c[k]))

# slices are not part of the protocol
q.json_results_append("Coordinate slice", probe(lambda: c[1:3]))

# assignment
c2 = q.Coordinate([1, 2, 3, 4])
c2[0] = 10
c2[-1] = 40
q.json_results_append("Coordinate setitem", np.array([float(v) for v in c2.to_list()]))

for key, val in ((1.5, 7), (4, 7)):
    q.json_results_append(f"Coordinate setitem key={key!r}", probe(lambda k=key, v=val: c2.__setitem__(k, v)))

# CoordinateD behaves the same
q.json_results_append("CoordinateD access", np.array([float(d[0]), float(d[-1])]))
for key in (4, -5, 1.5):
    q.json_results_append(f"CoordinateD d[{key!r}]", probe(lambda k=key: d[k]))

# conversions are unaffected
q.json_results_append("to_list", np.array([float(v) for v in c.to_list()]))
q.json_results_append("to_tuple", np.array([float(v) for v in c.to_tuple()]))
q.json_results_append("to_numpy", c.to_numpy())
q.json_results_append("iter", np.array([float(v) for v in c]))

# __len__ makes Coordinate a sized sequence, so NumPy converts it numerically
# instead of wrapping it in a 0-d object array
a_np = np.array(c)
a_as = np.asarray(c)
q.json_results_append("len(Coordinate)", float(len(c)))
q.json_results_append("len(CoordinateD)", float(len(d)))
q.json_results_append("np.array(Coordinate) dtype is not object", float(a_np.dtype != np.dtype(object)))
q.json_results_append("np.array(Coordinate)", np.array([float(v) for v in a_np]))
q.json_results_append("np.asarray(Coordinate) dtype is not object", float(a_as.dtype != np.dtype(object)))
q.json_results_append("np.asarray(Coordinate)", np.array([float(v) for v in a_as]))
q.json_results_append("np.asarray(CoordinateD)", np.array([float(v) for v in np.asarray(d)]))

# __bool__ is always True, so ``if c:`` never means "non-zero coordinate"
q.json_results_append("bool(non-zero Coordinate)", float(bool(c)))
q.json_results_append("bool(Coordinate())", float(bool(q.Coordinate())))
q.json_results_append("bool(CoordinateD())", float(bool(q.CoordinateD())))
q.json_results_append("Coordinate() == Coordinate()", float(q.Coordinate() == q.Coordinate()))

# __hash__ is consistent with __eq__, so coordinates can be dict keys / set members
c_same = q.Coordinate([1, 2, 3, 4])
d_same = q.CoordinateD([1.5, 2.5, 3.5, 4.5])
q.json_results_append("hash(Coordinate) is int", float(isinstance(hash(c), int)))
q.json_results_append("equal Coordinates hash equal", float(hash(c) == hash(c_same)))
q.json_results_append("equal CoordinateD hash equal", float(hash(d) == hash(d_same)))
q.json_results_append("unequal Coordinates hash differently", float(hash(c) != hash(q.Coordinate([0, 0, 0, 0]))))
q.json_results_append("Coordinate as set member", float(len({c, c_same, q.Coordinate([0, 0, 0, 0])}) == 2))
q.json_results_append("CoordinateD as set member", float(len({d, d_same, q.CoordinateD()}) == 2))
c_dict = {c: 1.0, q.Coordinate([0, 0, 0, 0]): 2.0}
q.json_results_append("Coordinate as dict key", np.array([c_dict[c_same], c_dict[q.Coordinate([0, 0, 0, 0])]]))
d_dict = {d: 1.0}
q.json_results_append("CoordinateD as dict key", float(d_dict[d_same]))

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-14)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
