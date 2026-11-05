#!/usr/bin/env python3

# Conformance test for the unified smear API.
#
# Asserts the public contract:
#   * every smear takes the iteration parameter `num_step` (last positional,
#     default 1) and returns a NEW object without touching its input,
#   * `prop_smear` / `prop_spatial_smear` keep `mom` (and the mode options)
#     after `num_step` so the historical positional order stays valid,
#   * `field_smear` spells the spatial switch `spatial_only` (keyword only),
#   * the pre-rename names are no longer exported, while the module names
#     `q.smear` and `q.smear_prop` are unchanged,
#   * `prop_spatial_smear_no_comm` is the only in-place (None-returning) smear,
#   * `field_smear(field, r, 2)` equals two applications of
#     `field_smear(field, r, 1)`.

import inspect

import numpy as np

import qlat as q

OLD_NAMES = [
    "gf_spatial_ape_smear",
    "gf_simple_no_proj_smear",
    "gf_block_stout_smear",
    "gf_local_stout_smear",
    "smear_field",
    "smear_field_step",
    "smear_field_step_local",
    "mk_smear_mom_kernel",
    "mk_spatial_smear_mom_kernel",
]

# Public smear operations: they all take `num_step` and return a new object.
NUM_STEP_NAMES = [
    "gf_ape_smear",
    "gf_ape_smear_spatial",
    "gf_ape_smear_no_proj",
    "gf_hyp_smear",
    "gf_stout_smear",
    "gf_stout_smear_block",
    "gf_stout_smear_local",
    "prop_smear",
    "prop_spatial_smear",
    "field_smear",
    "field_smear_step",
    "field_smear_step_local",
]

# These keep `mom` (and the mode options) after `num_step`.
TRAILING_OPTIONAL_NAMES = [
    "prop_smear",
    "prop_spatial_smear",
]

KERNEL_NAMES = [
    "mk_field_smear_mom_kernel",
    "mk_field_smear_mom_kernel_spatial",
]

def get_positional_params(func):
    sig = inspect.signature(func)
    return [
        param
        for param in sig.parameters.values()
        if param.kind
        in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        )
    ]

def check_num_step_signature(name, func, is_trailing_optional):
    positional = get_positional_params(func)
    names = [param.name for param in positional]
    assert "num_step" in names, f"{name}: no `num_step` parameter (got {names})"
    param = inspect.signature(func).parameters["num_step"]
    assert param.default == 1, (
        f"{name}: num_step default is {param.default!r}, expected 1"
    )
    if is_trailing_optional:
        index = names.index("num_step")
        assert names[index - 1] == "coef", (
            f"{name}: `num_step` does not follow `coef` (got {names})"
        )
        for extra in positional[index + 1 :]:
            assert extra.default is not inspect.Parameter.empty, (
                f"{name}: `{extra.name}` after num_step has no default"
            )
    else:
        assert positional[-1].name == "num_step", (
            f"{name}: last positional is `{positional[-1].name}`, expected `num_step`"
        )
    q.json_results_append(
        f"signature {name}: `num_step` iteration parameter, default 1"
    )

def get_field_sig(field):
    return np.array(q.get_data_sig(field, q.RngState()))

def get_local_field_sig(field):
    # Local-only signature: safe for the explicitly constructed single-node
    # geometry, whose fields are not part of the global node decomposition.
    return np.array(q.get_data_sig(np.asarray(field), q.RngState()))

def get_arr_sig(field):
    return np.array(q.get_data_sig_arr(field, q.RngState(), 3))

def check_returns_new_object(label, call, input_obj, sig_func):
    sig_before = sig_func(input_obj)
    result = call()
    sig_after = sig_func(input_obj)
    is_new = result is not None
    is_preserved = bool(np.array_equal(sig_before, sig_after))
    q.json_results_append(f"{label}: returns a new object = {is_new}")
    q.json_results_append(f"{label}: input field left unchanged = {is_preserved}")
    assert is_new, f"{label}: returned None, expected a new object"
    assert is_preserved, f"{label}: modified its input field"
    return result

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4],
]

q.begin_with_mpi(size_node_list)

rs = q.RngState("seed")
total_site = q.Coordinate([4, 4, 4, 4])
geo = q.Geometry(total_site)
# Rank-local geometry: `field_smear_step_local` works without communication
# whenever every site is rank-local, on any number of MPI ranks.
geo_local = q.Geometry(0, q.Coordinate([1, 1, 1, 1]), q.Coordinate([4, 4, 4, 4]))

# ---- signature contract ----

for name in NUM_STEP_NAMES:
    check_num_step_signature(name, getattr(q, name), name in TRAILING_OPTIONAL_NAMES)

prop_sig = inspect.signature(q.prop_smear)
q.json_results_append(
    f"signature prop_smear: mode default = {prop_sig.parameters['mode'].default!r}"
)
assert prop_sig.parameters["mode"].default == "qlat"

field_sig = inspect.signature(q.field_smear)
assert "spatial_only" in field_sig.parameters
assert "is_only_spatial" not in field_sig.parameters
q.json_results_append("signature field_smear: keyword `spatial_only` = True")

# ---- old names gone, module names kept ----

for name in OLD_NAMES:
    assert not hasattr(q, name), f"{name} is still exported by qlat"
q.json_results_append(f"old smear names absent = True ({len(OLD_NAMES)} checked)")

assert hasattr(q, "smear") and hasattr(q, "smear_prop")
q.json_results_append("module names `q.smear` and `q.smear_prop` kept = True")

# ---- gauge-field smears ----

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf"), 0.5, 10)

check_returns_new_object(
    "gf_ape_smear", lambda: q.gf_ape_smear(gf, 0.5, num_step=1), gf, get_arr_sig
)
check_returns_new_object(
    "gf_ape_smear_spatial",
    lambda: q.gf_ape_smear_spatial(gf, 0.5, num_step=1),
    gf,
    get_arr_sig,
)
check_returns_new_object(
    "gf_ape_smear_no_proj",
    lambda: q.gf_ape_smear_no_proj(gf, 0.5, num_step=1),
    gf,
    get_arr_sig,
)
check_returns_new_object(
    "gf_hyp_smear",
    lambda: q.gf_hyp_smear(gf, 0.75, 0.6, 0.3, num_step=1),
    gf,
    get_arr_sig,
)
check_returns_new_object(
    "gf_stout_smear", lambda: q.gf_stout_smear(gf, 0.1, num_step=1), gf, get_arr_sig
)
check_returns_new_object(
    "gf_stout_smear_block",
    lambda: q.gf_stout_smear_block(gf, q.Coordinate(), 0.1, num_step=1),
    gf,
    get_arr_sig,
)
check_returns_new_object(
    "gf_stout_smear_local",
    lambda: q.gf_stout_smear_local(gf, q.Coordinate(), 0.1, num_step=1),
    gf,
    get_arr_sig,
)

# ---- propagator smears ----

gf_ape = q.gf_ape_smear_spatial(gf, 0.5, num_step=3)
gf1 = q.mk_left_expanded_field(gf_ape)
prop = q.Prop(geo)
prop.set_rand(rs.split("prop"))
coef = 0.9375

check_returns_new_object(
    "prop_smear",
    lambda: q.prop_smear(prop, gf1, coef, num_step=1),
    prop,
    get_arr_sig,
)
check_returns_new_object(
    "prop_spatial_smear",
    lambda: q.prop_spatial_smear(prop, gf_ape, coef, num_step=1),
    prop,
    get_arr_sig,
)

# ---- field smears ----

f_c = q.FieldComplexD(geo, 3)
f_c.set_rand(rs.split("f_c"))
check_returns_new_object(
    "field_smear", lambda: q.field_smear(f_c, 2.0, num_step=1), f_c, get_field_sig
)

f_r = q.FieldRealD(geo, 1)
f_r.set_rand(rs.split("f_r"))
check_returns_new_object(
    "field_smear_step",
    lambda: q.field_smear_step(f_r, 0.5, num_step=1),
    f_r,
    get_field_sig,
)

f_r_local = q.FieldRealD(geo_local, 1)
f_r_local.set_rand(rs.split("f_r_local"))
check_returns_new_object(
    "field_smear_step_local",
    lambda: q.field_smear_step_local(f_r_local, 0.5, num_step=1),
    f_r_local,
    get_local_field_sig,
)

# ---- momentum-space kernels ----

total_site_tuple = total_site.to_tuple()
for name in KERNEL_NAMES:
    kernel = getattr(q, name)(total_site_tuple, 2.0)
    is_field = kernel is not None
    q.json_results_append(f"{name}: returns field = {is_field}")
    assert is_field, f"{name}: returned None"

# ---- the only in-place exception ----

geo_slice = q.Geometry(
    0, q.Coordinate([1, 1, 1, total_site[3]]), q.Coordinate([4, 4, 4, 1])
)
gf_local = q.GaugeField(geo_slice)
gf_local.set_rand(rs.split("gf_local"), 0.5, 10)
ff_list = [q.FermionField4d(geo_slice) for _ in range(2)]
for ff in ff_list:
    ff.set_rand(rs.split("ff"))
sig_before_list = [get_local_field_sig(ff) for ff in ff_list]
res_no_comm = q.prop_spatial_smear_no_comm(ff_list, gf_local, coef, 1)
sig_after_list = [get_local_field_sig(ff) for ff in ff_list]
is_none = res_no_comm is None
is_mutated = any(
    not np.array_equal(before, after)
    for before, after in zip(sig_before_list, sig_after_list)
)
q.json_results_append(f"prop_spatial_smear_no_comm: returns None = {is_none}")
q.json_results_append(f"prop_spatial_smear_no_comm: smears in place = {is_mutated}")
assert is_none and is_mutated
# Every other smear returned a new object above, so this one is the only
# in-place, None-returning exception.
q.json_results_append(
    "prop_spatial_smear_no_comm is the only in-place None-returning smear = True"
)

# ---- field_smear num_step composition ----

f_comp = q.FieldComplexD(geo, 1)
f_comp.set_rand(rs.split("f_comp"))
f_two = q.field_smear(f_comp, 2.0, num_step=2)
f_one_one = q.field_smear(q.field_smear(f_comp, 2.0, num_step=1), 2.0, num_step=1)
f_diff = f_two.copy()
f_diff -= f_one_one
diff_norm = float(f_diff.qnorm())
q.json_results_append(
    "field_smear num_step=2 vs two num_step=1 qnorm diff", diff_norm, 1e-10
)
q.json_results_append(f"field_smear num_step composition ok = {diff_norm < 1e-10}")
assert diff_norm < 1e-10, diff_norm

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-10)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
