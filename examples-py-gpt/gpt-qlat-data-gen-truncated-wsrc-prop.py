#!/usr/bin/env python3

"""
Generate light and strange quark wall-source propagators with time-direction
truncated gauge fields.
Detailed Description
====================
This script generates wall-source propagators for light and strange quarks
using gauge fields that are truncated in the time direction. For each source
time slice, a sub-volume of the gauge field is extracted centered on that
source. The Dirac inversion is performed on this truncated gauge field, and
only propagator data within the truncated region is saved.
The key idea is that instead of inverting on the full lattice gauge field,
we restrict the gauge field to a window of size ``2 * t_half + 1`` time
slices centered on the wall-source time slice, padded to a multiple of
``field_trunc_t_size_divisor``. This can be useful for:
- Studying finite-volume effects in the time direction
- Reducing computational cost for localized measurements
- Isolating propagation within a specific time region
For each source time slice ``t_src``, the truncated gauge field covers the
time range ``[t_src - t_half, t_src + t_half + 1]`` (modulo ``t_size``). The
source is placed at time index ``t_half`` in the truncated geometry.
One type of propagator sink is saved:
1. **Wall sink** (``psel-prop-wsrc-trunc-{flavor}``): The propagator summed over
   all spatial sites at each time slice within the truncated region. This is
   stored as a ``PselProp`` with one entry per time slice of the truncated
   geometry, saved in ``.lat`` format.
The output data layout is::
    {job_tag}/psel-prop-wsrc-trunc-{flavor}/traj-{traj}/
        wsnk.lat      # Wall-sink PselProp data per (idx, tslice, inv_type, inv_acc)
Parameters (via ``set_param``)
------------------------------\n
- ``measurement.field_trunc_half_width_list``: (Mandatory) List of half-widths
  of the truncated time region. The calculation is repeated for each
  half-width. The full truncated region spans ``2 * half_width + 1`` time
  slices, padded to a multiple of ``field_trunc_t_size_divisor``.
- ``measurement.field_trunc_t_size_divisor``: (Mandatory) The truncated time
  extent must be a multiple of this value. Padding slices are added beyond
  ``2 * half_width + 1`` as needed.\n
Test Mode
---------\n
With ``--test`` flag, a small ``test-4nt16-checker`` configuration is used
(lattice 4x4x4x16, 2 MPI ranks along time). The test generates a sample
gauge field via HMC (5 trajectories, 8 Wilson flow steps), then computes
truncated wall-source propagators for light and strange quarks with
``field_trunc_half_width_list = [3]``, ``field_trunc_t_size_divisor = 8``
(producing a truncated geometry of 4x4x4x8).
All time slices are used as source time slices.\n
Usage
-----\n
::\n
    mpirun -n <N> python3 gpt-qlat-data-gen-truncated-wsrc-prop.py
    mpirun -n <N> python3 gpt-qlat-data-gen-truncated-wsrc-prop.py --job_tag_list 24D,32Dfine\n
See Also
--------\n
- ``gpt-qlat-data-gen-eta-c.py``: Similar wall-source propagator generation
  for charm quarks (without truncation)
- ``qlat_scripts.v1.gen_data``: Core propagator generation utilities
"""

import qlat_gpt as qg
import qlat as q
import numpy as np

import qlat_scripts.v1 as qs

from qlat_scripts.v1 import (
    load_path_list,
    get_param,
    set_param,
    run_params,
    check_job,
    run_gf,
    run_gt,
    get_save_path,
    is_test,
)

### ------

load_path_list[:] = [
    "results",
    "qcddata",
    "qcddata-1",
    "qcddata-2",
    "qcddata-3",
    "qcddata-4",
    "qcddata-5",
    "qcddata-6",
    "qcddata-7",
    "qcddata-8",
    "qcddata-9",
]

### ------

@q.timer
def mk_field_truncated(field, t_start, t_end):
    """
    Create a field truncated in the time direction.\n
    Extracts a sub-volume of the field covering the time range
    ``[t_start, t_end)`` (modulo ``t_size``). The truncated geometry has the
    same spatial extent as the original but with time extent
    ``t_size_trunc`` (rounded up to be divisible by the number of
    time-direction MPI nodes).\n
    The field data is copied from the original field (and its time-shifted
    copies) for all sites whose global time coordinate falls within the
    truncated range. Sites that are not local in the truncated geometry
    are skipped.\n
    Parameters
    ----------
    field : q.FieldBase
        The input field (e.g. ``GaugeField``, ``Prop``, ``FieldComplexD``).
    t_start : int
        Start of the truncated time range (inclusive, in the original
        geometry's global time coordinates).
    t_end : int
        End of the truncated time range (exclusive, in the original
        geometry's global time coordinates). Wraps around modulo ``t_size``.\n
    Returns
    -------
    field_trunc : same type as ``field``
        Field on the truncated geometry with data copied from the original.
    """
    geo = field.geo
    total_site = geo.total_site
    t_size = total_site[3]
    t_size_trunc = (t_end - t_start) % t_size
    if t_size_trunc == 0:
        t_size_trunc = t_size
    assert t_size_trunc <= t_size
    t_size_trunc_valid = t_size_trunc
    size_node_t = geo.size_node[3]
    t_size_trunc = ((t_size_trunc + size_node_t - 1) // size_node_t) * size_node_t
    total_site_trunc = q.Coordinate(
        [total_site[0], total_site[1], total_site[2], t_size_trunc]
    )
    geo_trunc = q.Geometry(total_site_trunc)
    field_trunc = type(field)(geo_trunc, field.multiplicity)
    t_offset = t_start % t_size
    node_site_t = geo.node_site[3]
    n_shift = geo.size_node[3]
    shift = q.Coordinate([0, 0, 0, node_site_t])
    node_site_trunc = geo_trunc.node_site
    coor_node = geo.coor_node
    current_field = field.copy()
    current_xg_field = q.mk_xg_field(geo)
    field_trunc_arr = np.asarray(field_trunc)
    for i_shift in range(n_shift):
        current_field_arr = np.asarray(current_field)
        xg_arr = current_xg_field[:]
        for index in range(current_field.geo.local_volume):
            xg = xg_arr[index]
            t_orig = xg[3] % t_size
            t_local = (t_orig - t_offset) % t_size
            if t_local >= t_size_trunc_valid:
                continue
            xl = [xg[i] - coor_node[i] * node_site_trunc[i] for i in range(3)]
            xl3 = t_local - coor_node[3] * node_site_trunc[3]
            if xl3 < 0 or xl3 >= node_site_trunc[3]:
                continue
            index_trunc = xl[0] + node_site_trunc[0] * (
                xl[1] + node_site_trunc[1] * (xl[2] + node_site_trunc[2] * xl3)
            )
            field_trunc_arr[index_trunc] = current_field_arr[index]
        if i_shift < n_shift - 1:
            current_field = current_field.shift(shift)
            current_xg_field = current_xg_field.shift(shift)
    return field_trunc

@q.timer
def mk_gf_truncated(gf, t_center, t_half, t_size_divisor=1):
    """
    Create a gauge field truncated in the time direction.
    Extracts a sub-volume of the gauge field centered at ``t_center`` with
    half-width ``t_half``. The truncated geometry has the same spatial extent
    as the original but with time extent ``2 * t_half + 1``, padded to be a
    multiple of ``t_size_divisor``.
    Implemented via :func:`mk_field_truncated` with
    ``t_start = (t_center - t_half) % t_size`` and
    ``t_end = (t_start + t_size_trunc) % t_size``.
    Parameters
    ----------
    gf : q.GaugeField
        The full gauge field.
    t_center : int
        Center time slice of the truncated region (in the original geometry).
    t_half : int
        Half-width of the truncated region. The needed truncated time extent
        is ``2 * t_half + 1``, padded to a multiple of
        ``t_size_divisor``.
    t_size_divisor : int
        The truncated time extent must be a multiple of this value. Padding
        slices are added beyond ``2 * t_half + 1`` as needed.
    Returns
    -------
    gf_trunc : q.GaugeField
        Gauge field on the truncated geometry.
    t_offset : int
        The global time index corresponding to time index 0 of the truncated
        geometry. ``t_offset = (t_center - t_half) % t_size``.
    t_size_trunc : int
        The actual truncated time extent (padded to a multiple of
        ``t_size_divisor``).
    """
    total_site = gf.geo.total_site
    t_size = total_site[3]
    t_size_trunc = 2 * t_half + 1
    r = t_size_trunc % t_size_divisor
    if r != 0:
        t_size_trunc += t_size_divisor - r
    assert t_size_trunc <= t_size
    t_start = (t_center - t_half) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gf_trunc = mk_field_truncated(gf, t_start, t_end)
    geo_trunc = gf_trunc.geo
    tslice_target_list = list(range(2 * t_half, t_size_trunc))
    gf_arr = np.asarray(gf_trunc)
    xg_arr = q.mk_xg_field(geo_trunc)[:]
    for index in range(geo_trunc.local_volume):
        xg = xg_arr[index]
        xg_t = xg[3]
        if xg_t in tslice_target_list:
            gf_arr[index, 3, :, :] = 0
    return gf_trunc, t_start, t_size_trunc

@q.timer
def mk_gt_truncated(gt, t_center, t_half, t_size_divisor=1):
    """
    Create a gauge transform field truncated in the time direction.
    Extracts a sub-volume of the gauge transform centered at ``t_center``
    with half-width ``t_half``. The truncated geometry has the same spatial
    extent as the original but with time extent ``2 * t_half + 1``, padded to
    be a multiple of ``t_size_divisor``.
    Parameters
    ----------
    gt : q.GaugeTransform
        The full gauge transform field.
    t_center : int
        Center time slice of the truncated region (in the original geometry).
    t_half : int
        Half-width of the truncated region.
    t_size_divisor : int
        The truncated time extent must be a multiple of this value.
    Returns
    -------
    gt_trunc : q.GaugeTransform
        Gauge transform on the truncated geometry.
    t_offset : int
        The global time index corresponding to time index 0 of the truncated
        geometry. ``t_offset = (t_center - t_half) % t_size``.
    t_size_trunc : int
        The actual truncated time extent (padded to a multiple of
        ``t_size_divisor``).
    """
    total_site = gt.geo.total_site
    t_size = total_site[3]
    t_size_trunc = 2 * t_half + 1
    r = t_size_trunc % t_size_divisor
    if r != 0:
        t_size_trunc += t_size_divisor - r
    assert t_size_trunc <= t_size
    t_start = (t_center - t_half) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gt_trunc = mk_field_truncated(gt, t_start, t_end)
    return gt_trunc, t_start, t_size_trunc

@q.timer
def mk_selected_points_truncated(sp, idx_start, idx_end):
    """
    Return a copy of ``sp[idx_start:idx_end]`` as a new SelectedPoints.
    """
    n_keep = idx_end - idx_start
    total_site = sp.psel.total_site
    psel_sub = q.PointsSelection(total_site, n_keep)
    np.asarray(psel_sub)[:] = np.asarray(sp.psel)[idx_start:idx_end]
    sp_trunc = type(sp)(psel_sub)
    sp_trunc @= sp
    return sp_trunc

### ------

@q.timer(is_timer_fork=True)
def run_prop_wsrc_truncated_save(job_tag, traj, *, get_gf, get_gt, inv_type):
    """
    Compute and save truncated wall-source propagators for one quark flavor.\n
    Saves wall-sink data as psel in .lat format.
    """
    fname = q.get_fname()
    inv_type_name = ["light", "strange"][inv_type]
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_site[3]
    t_half_list = get_param(job_tag, "measurement", "field_trunc_half_width_list")
    t_size_divisor = get_param(job_tag, "measurement", "field_trunc_t_size_divisor")
    tslice_list = get_truncated_wsrc_tslice_list(job_tag, traj)
    gf = get_gf()
    gt = get_gt()
    inv_acc = get_param(job_tag, "measurement", "field_trunc_inv_acc_list")[inv_type]
    acc_tag = f"accuracy={inv_acc}"
    for t_half in t_half_list:
        path_ws = f"{job_tag}/psel-prop-wsrc-trunc-{inv_type_name}/traj-{traj}/t_half-{t_half}"
        if q.does_file_exist_qar_sync_node(get_save_path(path_ws + "/checkpoint.txt")):
            continue
        if not q.obtain_lock(
            f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}-t_half-{t_half}"
        ):
            continue
        qar_ws = q.open_qar_info(get_save_path(path_ws + ".qar"), "a")
        for idx, tslice in enumerate(tslice_list):
            tag_base = f"tslice={tslice} ; t_half={t_half} ; type={inv_type}"
            tag = f"{tag_base} ; {acc_tag}"
            if q.bcast_py(qar_ws.has_regular_file(f"idx={idx} ; {tag} ; wsnk.lat")):
                q.displayln_info(f"WARNING: qar_ws has: old tag {tag}")
                continue
            has_file = q.bcast_py(qar_ws.has_regular_file(f"{tag} ; wsnk.lat"))
            if has_file:
                continue
            q.check_time_limit()
            gf_trunc, t_offset, t_size_trunc = mk_gf_truncated(
                gf, tslice, t_half, t_size_divisor
            )
            gt_trunc, t_offset_gt, t_size_trunc_gt = mk_gt_truncated(
                gt, tslice, t_half, t_size_divisor
            )
            assert t_offset == t_offset_gt
            assert t_size_trunc == t_size_trunc_gt
            geo_trunc = gf_trunc.geo
            t_src_trunc = t_half
            inv_trunc = qs.get_inv(
                gf_trunc, job_tag, inv_type, inv_acc, gt=gt_trunc, eig=None
            )
            q.displayln_info(0, f"truncated gf total_site: {gf_trunc.geo.total_site}")
            src = q.mk_wall_src(geo_trunc, t_src_trunc)
            sol_trunc = inv_trunc * src
            ps_prop_ws = sol_trunc.glb_sum_tslice()
            ps_prop_ws = mk_selected_points_truncated(ps_prop_ws, 0, 2 * t_half + 1)
            q.json_results_append(f"sol_trunc qnorm {tag}", sol_trunc.qnorm(), 1e-6)
            q.json_results_append(
                f"sol_trunc sig {tag}", q.get_data_sig(sol_trunc, q.RngState())
            )
            q.json_results_append(f"ps_prop_ws qnorm {tag}", ps_prop_ws.qnorm(), 1e-6)
            q.json_results_append(
                f"ps_prop_ws sig {tag}", q.get_data_sig(ps_prop_ws, q.RngState())
            )
            qar_ws.write(
                f"{tag} ; wsnk.lat", "", ps_prop_ws.save_str(), skip_if_exist=True
            )
            qar_ws.flush()
            q.clean_cache(q.cache_inv)
        qar_ws.write("checkpoint.txt", "", "", skip_if_exist=True)
        qar_ws.flush()
        qar_ws.close()
        q.release_lock()

### ------

@q.cache_call()
@q.timer()
def get_truncated_wsrc_tslice_list(job_tag, traj):
    """
    Generate the list of source time slices for truncated wall-source
    propagators.\n
    Returns all time slices ``[0, 1, ..., t_size-1]``.
    """
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    return list(range(t_size))

### ------

@q.timer(is_timer_fork=True)
def run_job(job_tag, traj):
    traj_gf = traj
    if job_tag[:5] == "test-":
        traj_gf = 1000
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_site[3]
    t_half_list = get_param(job_tag, "measurement", "field_trunc_half_width_list")
    fns_produce = []
    for inv_type_name in ["light", "strange"]:
        for t_half in t_half_list:
            fns_produce.append(
                f"{job_tag}/psel-prop-wsrc-trunc-{inv_type_name}/traj-{traj}/t_half-{t_half}/checkpoint.txt"
            )
    fns_need = [
        f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
        (
            f"{job_tag}/configs/ckpoint_lat.{traj_gf}",
            f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",
        ),
        f"{job_tag}/psel-prop-psrc-strange/traj-{traj}/checkpoint.txt",
    ]
    if job_tag[:5] == "test-":
        fns_need = []
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    for inv_type in [0, 1]:
        run_prop_wsrc_truncated_save(
            job_tag,
            traj,
            get_gf=get_gf,
            get_gt=get_gt,
            inv_type=inv_type,
        )

### ------

job_tag = "24D"
set_param(job_tag, "traj_list")(list(range(500, 5000, 10)))
set_param(job_tag, "measurement", "field_trunc_half_width_list")([1, 3, 5, 7])
set_param(job_tag, "measurement", "field_trunc_t_size_divisor")(4)
set_param(job_tag, "measurement", "field_trunc_inv_acc_list")([2, 2])
set_param(job_tag, "cg_params-0-2", "maxiter")(300)
set_param(job_tag, "cg_params-0-2", "maxcycle")(3)
set_param(job_tag, "cg_params-1-2", "maxiter")(300)
set_param(job_tag, "cg_params-1-2", "maxcycle")(1)

job_tag = "32Dfine"
set_param(job_tag, "traj_list")(list(range(500, 5000, 10)))
set_param(job_tag, "measurement", "field_trunc_half_width_list")([1, 3, 5, 7, 9])
set_param(job_tag, "measurement", "field_trunc_t_size_divisor")(4)
set_param(job_tag, "measurement", "field_trunc_inv_acc_list")([2, 2])
set_param(job_tag, "cg_params-0-2", "maxiter")(300)
set_param(job_tag, "cg_params-0-2", "maxcycle")(3)
set_param(job_tag, "cg_params-1-2", "maxiter")(300)
set_param(job_tag, "cg_params-1-2", "maxcycle")(1)

job_tag = "64I"
set_param(job_tag, "traj_list")(list(range(1200, 3680, 80)))
set_param(job_tag, "measurement", "field_trunc_half_width_list")([1, 5, 9, 13, 17])
set_param(job_tag, "measurement", "field_trunc_t_size_divisor")(4)
set_param(job_tag, "measurement", "field_trunc_inv_acc_list")([2, 2])
set_param(job_tag, "cg_params-0-2", "maxiter")(300)
set_param(job_tag, "cg_params-0-2", "maxcycle")(3)
set_param(job_tag, "cg_params-1-2", "maxiter")(300)
set_param(job_tag, "cg_params-1-2", "maxcycle")(1)

# ----

job_tag = "test-4nt16-checker"
#
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
#
set_param(job_tag, "seed")("test-4nt16-checker")
set_param(job_tag, "total_site")(
    [
        4,
        4,
        4,
        16,
    ]
)
set_param(job_tag, "load_config_params", "twist_boundary_at_boundary")(
    [
        0.0,
        0.0,
        0.0,
        -0.5,
    ]
)
#
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(2)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
#
set_param(job_tag, "quark_flavor_list")(
    [
        "light",
        "strange",
    ]
)
set_param(job_tag, "quark_mass_list")(
    [
        0.01,
        0.04,
    ]
)
set_param(job_tag, "fermion_params", 0, 0)(
    {
        "Ls": 8,
        "M5": 1.8,
        "b": 1.5,
        "c": 0.5,
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
    }
)
for inv_type, mass in enumerate(get_param(job_tag, "quark_mass_list")):
    set_param(job_tag, "fermion_params", inv_type, 0)(
        get_param(job_tag, "fermion_params", 0, 0).copy()
    )
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [0, 1, 2]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(
            get_param(job_tag, "fermion_params", inv_type, 0).copy()
        )
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(10)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(1 + inv_acc)
#
set_param(job_tag, "measurement", "field_trunc_half_width_list")([1, 2, 3])
set_param(job_tag, "measurement", "field_trunc_t_size_divisor")(8)
set_param(job_tag, "measurement", "field_trunc_inv_acc_list")([2, 2])
set_param(job_tag, "cg_params-0-2", "maxcycle")(1)
set_param(job_tag, "cg_params-1-2", "maxcycle")(1)

# ----

##################### CMD options #####################

job_tag_list_default = [
    "test-4nt16-checker",
]
job_tag_list_str_default = ",".join(job_tag_list_default)
job_tag_list = q.get_arg("--job_tag_list", default=job_tag_list_str_default).split(",")

random_order = q.get_arg("--random-order", default=None) is not None

#######################################################

def gracefully_finish():
    q.displayln_info("Begin to gracefully_finish.")
    q.timer_display()
    if is_test():
        q.json_results_append(
            f"q.obtained_lock_history_list={q.obtained_lock_history_list}"
        )
        q.check_log_json(__file__)
    qg.end_with_gpt()
    q.displayln_info("CHECK: finished successfully.")
    exit()

def try_gracefully_finish():
    """
    Call `gracefully_finish` if not test and if some work is done
    (q.obtained_lock_history_list != []).
    """
    if (not is_test()) and (len(q.obtained_lock_history_list) > 0):
        gracefully_finish()

if __name__ == "__main__":
    qg.begin_with_gpt()
    q.check_time_limit()
    #
    job_tag_traj_list = []
    for job_tag in job_tag_list:
        run_params(job_tag)
        traj_list = get_param(job_tag, "traj_list")
        for traj in traj_list:
            job_tag_traj_list.append(
                (
                    job_tag,
                    traj,
                )
            )
    if not is_test() and random_order:
        job_tag_traj_list = q.random_permute(
            job_tag_traj_list, q.RngState(f"{q.get_time()}")
        )
        job_tag_traj_list = q.get_comm().bcast(job_tag_traj_list)
    for job_tag, traj in job_tag_traj_list:
        q.check_time_limit()
        run_job(job_tag, traj)
        q.clean_cache()
        try_gracefully_finish()
    #
    gracefully_finish()

# ----
