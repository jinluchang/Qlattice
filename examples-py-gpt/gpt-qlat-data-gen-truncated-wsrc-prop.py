#!/usr/bin/env python3

"""
Generate light and strange quark wall-source propagators with time-direction
truncated gauge fields.\n
Detailed Description
====================\n
This script generates wall-source propagators for light and strange quarks
using gauge fields that are truncated in the time direction. For each source
time slice, a sub-volume of the gauge field is extracted centered on that
source. The Dirac inversion is performed on this truncated gauge field, and
only propagator data within the truncated region is saved.\n
The key idea is that instead of inverting on the full lattice gauge field,
we restrict the gauge field to a window of size ``2 * t_half + 1`` time
slices centered on the wall-source time slice. This can be useful for:\n
- Studying finite-volume effects in the time direction
- Reducing computational cost for localized measurements
- Isolating propagation within a specific time region\n
For each source time slice ``t_src``, the truncated gauge field covers the
time range ``[t_src - t_half, t_src + t_half]`` (modulo ``t_size``). The
source is placed at the center of this truncated region (time index
``t_half`` in the truncated geometry).\n
Two types of propagator sinks are saved:\n
1. **Wall sink** (``psel-prop-wsrc-trunc-{flavor}``): The propagator summed over
   all spatial sites at each time slice within the truncated region. This is
   stored as a ``PselProp`` with one entry per time slice of the truncated
   geometry, saved in ``.lat`` format.\n
The output data layout is::\n
    {job_tag}/psel-prop-wsrc-trunc-{flavor}/traj-{traj}/
        wsnk.lat      # Wall-sink PselProp data per (idx, tslice, inv_type, inv_acc)\n
Parameters (via ``set_param``)
------------------------------\n
- ``measurement.num_wsrc_trunc``: Number of source time slices per trajectory
- ``measurement.wsrc_trunc_half_width``: Half-width of the truncated time
  region (default ``t_size // 4``). The full truncated region spans
  ``2 * half_width + 1`` time slices.
- ``measurement.wsrc_trunc_acc_list``: List of inversion accuracy levels to
  compute (default ``[0, 2]`` for sloppy and exact)\n
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

import qlat_scripts.v1 as qs

from qlat_scripts.v1 import (
    load_path_list,
    get_param,
    set_param,
    get_job_seed,
    run_params,
    check_job,
    run_gf,
    run_gt,
    get_load_path,
    get_save_path,
    is_test,
)

### ------

load_path_list[:] = [
    "results",
    "qcddata",
    "/lustre/orion/lgt119/proj-shared/ljin/qcddata4",
    "/lustre/orion/lgt119/proj-shared/ljin/qcddata5",
    "/lustre20/volatile/qcdqedta/qcddata",
    "/lustre20/volatile/decay0n2b/qcddata",
    "/lustre20/volatile/pqpdf/ljin/qcddata",
    "/lustre20/volatile/decay0n2b/qcddata/qcddata4",
    "/lustre20/volatile/decay0n2b/qcddata/qcddata3",
    "/lustre20/volatile/decay0n2b/qcddata/qcddata1",
    "/data1/qcddata4",
    "/data1/qcddata3",
    "/data2/qcddata3-prop",
    "/data1/qcddata1",
]

### ------

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
    current_field = field
    current_xg_field = q.mk_xg_field(geo)
    for i_shift in range(n_shift):
        xg_arr = current_xg_field[:]
        for index in range(current_field.geo.local_volume):
            xg = xg_arr[index]
            t_local = (xg[3] - t_offset) % t_size
            if t_local >= t_size_trunc:
                continue
            xl = [xg[i] - coor_node[i] * node_site_trunc[i] for i in range(3)]
            xl3 = t_local - coor_node[3] * node_site_trunc[3]
            if xl3 < 0 or xl3 >= node_site_trunc[3]:
                continue
            index_trunc = xl[0] + node_site_trunc[0] * (
                xl[1] + node_site_trunc[1] * (xl[2] + node_site_trunc[2] * xl3)
            )
            field_trunc.get_elems(index_trunc)[:] = current_field.get_elems(index)
        if i_shift < n_shift - 1:
            current_field = q.field_shift(current_field, shift)
            current_xg_field = q.field_shift(current_xg_field, shift)
    return field_trunc

def mk_gf_truncated(gf, t_center, t_half):
    """
    Create a gauge field truncated in the time direction.\n
    Extracts a sub-volume of the gauge field centered at ``t_center`` with
    half-width ``t_half``. The truncated geometry has the same spatial extent
    as the original but with time extent ``2 * t_half + 1``.\n
    Implemented via :func:`mk_field_truncated` with
    ``t_start = (t_center - t_half) % t_size`` and
    ``t_end = (t_center + t_half + 2) % t_size``.\n
    Parameters
    ----------
    gf : q.GaugeField
        The full gauge field.
    t_center : int
        Center time slice of the truncated region (in the original geometry).
    t_half : int
        Half-width of the truncated region. The truncated time extent is
        ``2 * t_half + 1``.\n
    Returns
    -------
    gf_trunc : q.GaugeField
        Gauge field on the truncated geometry.
    t_offset : int
        The global time index corresponding to time index 0 of the truncated
        geometry. ``t_offset = (t_center - t_half) % t_size``.
    """
    total_site = gf.geo.total_site
    t_size = total_site[3]
    t_start = (t_center - t_half) % t_size
    t_end = (t_center + t_half + 2) % t_size
    gf_trunc = mk_field_truncated(gf, t_start, t_end)
    return gf_trunc, t_start

### ------

@q.timer_verbose
def run_truncated_wsrc_prop(job_tag, traj, get_gf, get_gt):
    """
    Main driver: generate truncated wall-source propagators for light and
    strange quarks and save wall-sink data.
    """
    fname = q.get_fname()
    fn_checkpoint = f"{job_tag}/truncated-wsrc-prop/traj-{traj}/checkpoint.txt"
    if get_load_path(fn_checkpoint) is not None:
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    t_half = get_param(
        job_tag, "measurement", "wsrc_trunc_half_width", default=t_size // 4
    )
    tslice_list = get_truncated_wsrc_tslice_list(job_tag, traj)
    for inv_type in [0, 1]:
        run_prop_wsrc_truncated_save(
            job_tag,
            traj,
            get_gf=get_gf,
            get_gt=get_gt,
            inv_type=inv_type,
            tslice_list=tslice_list,
            t_half=t_half,
        )
    q.qtouch_info(get_save_path(fn_checkpoint))
    q.release_lock()

@q.timer_verbose
def run_prop_wsrc_truncated_save(
    job_tag, traj, *, get_gf, get_gt, inv_type, tslice_list, t_half
):
    """
    Compute and save truncated wall-source propagators for one quark flavor.\n
    Saves wall-sink data as psel in .lat format.
    """
    inv_type_name = ["light", "strange"][inv_type]
    gf = get_gf()
    inv_acc = 2
    acc_tag = f"accuracy={inv_acc}"
    path_ws = f"{job_tag}/psel-prop-wsrc-trunc-{inv_type_name}/traj-{traj}"
    qar_ws = q.open_qar_info(get_save_path(path_ws + ".qar"), "a")
    for idx, tslice in enumerate(tslice_list):
        tag_base = f"idx={idx} ; tslice={tslice} ; type={inv_type}"
        tag = f"{tag_base} ; {acc_tag}"
        if qar_ws.has_regular_file(f"{tag} ; wsnk.lat"):
            continue
        gf_trunc, t_offset = mk_gf_truncated(gf, tslice, t_half)
        geo_trunc = gf_trunc.geo
        t_src_trunc = t_half
        inv_trunc = qs.get_inv(gf_trunc, job_tag, inv_type, inv_acc, gt=None, eig=None)
        src = q.mk_wall_src(geo_trunc, t_src_trunc)
        sol_trunc = inv_trunc * src
        ps_prop_ws = sol_trunc.glb_sum_tslice()
        qar_ws.write(f"{tag} ; wsnk.lat", "", ps_prop_ws.save_str(), skip_if_exist=True)
        qar_ws.flush()
    qar_ws.write("checkpoint.txt", "", "", skip_if_exist=True)
    qar_ws.flush()
    qar_ws.close()
    q.clean_cache(q.cache_inv)

### ------

@q.cache_call()
@q.timer_verbose
def get_truncated_wsrc_tslice_list(job_tag, traj):
    """
    Generate the list of source time slices for truncated wall-source
    propagators.\n
    The time slices are approximately evenly spaced across the temporal
    extent of the lattice, with the starting point randomized per
    trajectory for unbiased sampling.
    """
    num_wsrc = get_param(job_tag, "measurement", "num_wsrc_trunc", default=2)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    t_size = total_site[3]
    rs = q.RngState(f"{get_job_seed(job_tag)}-{traj}-get_truncated_wsrc_tslice_list")
    t_start = rs.u_rand_gen() * t_size
    t_sep = t_size / num_wsrc
    tslice_list = [round(t_start + i * t_sep) % t_size for i in range(num_wsrc)]
    return tslice_list

### ------

@q.timer(is_timer_fork=True)
def run_job(job_tag, traj):
    traj_gf = traj
    if job_tag[:5] == "test-":
        traj_gf = 1000
    fns_produce = [
        f"{job_tag}/truncated-wsrc-prop/traj-{traj}/checkpoint.txt",
    ]
    fns_need = [
        f"{job_tag}/gauge-transform/traj-{traj_gf}.field",
        (
            f"{job_tag}/configs/ckpoint_lat.{traj_gf}",
            f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",
        ),
    ]
    if job_tag[:5] == "test-":
        fns_need = []
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    run_truncated_wsrc_prop(job_tag, traj, get_gf, get_gt)

### ------

job_tag = "24D"
set_param(job_tag, "traj_list")(list(range(1000, 5100, 80)))
set_param(job_tag, "measurement", "num_wsrc_trunc")(2)
set_param(job_tag, "measurement", "wsrc_trunc_half_width")(8)

job_tag = "32Dfine"
set_param(job_tag, "traj_list")(list(range(520, 2600, 40)))
set_param(job_tag, "measurement", "num_wsrc_trunc")(2)
set_param(job_tag, "measurement", "wsrc_trunc_half_width")(10)

job_tag = "64I"
set_param(job_tag, "traj_list")(list(range(1200, 3680, 80)))
set_param(job_tag, "measurement", "num_wsrc_trunc")(2)
set_param(job_tag, "measurement", "wsrc_trunc_half_width")(16)

# ----

job_tag = "test-4nt8-checker"
#
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
#
set_param(job_tag, "seed")("test-4nt8-checker")
set_param(job_tag, "total_site")(
    [
        4,
        4,
        4,
        8,
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
set_param(job_tag, "measurement", "num_wsrc_trunc")(2)
set_param(job_tag, "measurement", "wsrc_trunc_half_width")(2)

# ----

##################### CMD options #####################

job_tag_list_default = [
    "test-4nt8-checker",
]
job_tag_list_str_default = ",".join(job_tag_list_default)
job_tag_list = q.get_arg("--job_tag_list", default=job_tag_list_str_default).split(",")

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
    if not is_test():
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
