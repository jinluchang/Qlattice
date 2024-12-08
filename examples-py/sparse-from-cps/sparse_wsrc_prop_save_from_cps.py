import numpy as np
import qlat as q
import qlat_cps as qc
import qlat_gpt as qg
import gpt as g

from qlat_scripts.v1 import *
from sparse_prop_selection import *

@q.cache_call()
@q.timer
def get_fn_gauge_transform_avail_from_cps(job_tag, traj):
    """
    Check a Coulumb gauge transformation matrix is available.
    Return `fn` is the gauge transformation matrix is available.
    Otherwise, return `None`.
    #
    TODO: add actual CPS data location.
    """
    path = f"{job_tag}/gauge-transform-cps/traj-{traj}.gfix"
    fn = get_load_path(path)
    if fn:
        return fn
    #
    # TODO: add code here to return actual CPS data filename (full path) or `None`.
    #
    return False

@q.cache_call()
@q.timer_verbose
def check_prop_wsrc_avail_from_cps(
        job_tag, traj, *, inv_type=None,
    ):
    """
    Check the data is available.
    If `inv_type is None`, then return `True` if any `inv_type` data is available.
    #
    TODO: add actual CPS data check.
    """
    if inv_type is None:
        for inv_type in [ 0, 1, 2, ]:
            if check_prop_wsrc_avail_from_cps(job_tag, traj, inv_type=inv_type):
                return True
    inv_type_name_list = [ "light", "strange", "charm", ]
    inv_type_name = inv_type_name_list[inv_type]
    path = f"{job_tag}/prop-wsrc-full-cps-{inv_type_name}/traj-{traj}"
    if get_load_path(f"{path}/checkpoint.txt"):
        return True
    #
    # TODO: add actual CPS data check here.
    #
    return False

@q.cache_call()
@q.timer
def get_fn_prop_wsrc_avail_from_cps(
        job_tag, traj, *, tslice, inv_type, inv_acc,
    ):
    """
    Check a single propagator is available.
    Return `fn` is the propagator is available.
    Otherwise, return `None`.
    #
    TODO: add actual CPS data location.
    """
    inv_type_name_list = [ "light", "strange", "charm", ]
    inv_type_name = inv_type_name_list[inv_type]
    path = f"{job_tag}/prop-wsrc-full-cps-{inv_type_name}/traj-{traj}"
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    fn = get_load_path(f"{path}/{tag}.cps-prop")
    if fn:
        return fn
    #
    # TODO: add code here to return actual CPS data filename (full path) or `None`.
    #
    return None

@q.timer_verbose
def run_gauge_transform_from_cps(job_tag, traj):
    fname = q.get_fname()
    is_doing_something = None
    tfn = f"{job_tag}/gauge-transform/traj-{traj}.field"
    tfn_cps = f"{job_tag}/gauge-transform/traj-{traj}.gfix"
    path_gt = get_load_path(tfn)
    if path_gt is not None:
        assert get_load_path(tfn_cps) is not None
        q.displayln_info(f"{fname}: {job_tag}/{traj} already done.")
        return None
    if not check_prop_wsrc_avail_from_cps(job_tag, traj):
        q.displayln_info(f"{fname}: {job_tag}/{traj} CPS data not available.")
        return None
    q.check_stop()
    q.check_time_limit()
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    with q.TimerFork(max_call_times_for_always_show_info=0, verbose=1, show_display=True):
        fn = get_fn_gauge_transform_avail_from_cps(job_tag, traj)
        assert fn is not None # We need to have this data if we have propagators
        total_site = q.Coordinate(get_param(job_tag, "total_site"))
        geo = q.Geometry(total_site)
        gt = q.GaugeTransform(geo)
        gt.load_cps(fn)
        gt.save_cps(get_save_path(tfn_cps))
        gt.save_double(get_save_path(tfn))
        is_doing_something = True
        q.release_lock()
    return is_doing_something

@q.timer_verbose
def run_prop_wsrc_sparse_from_cps(job_tag, traj, *, inv_type, get_fsel, get_psel):
    fname = q.get_fname()
    is_doing_something = None
    inv_type_name_list = [ "light", "strange", "charm", ]
    inv_type_name = inv_type_name_list[inv_type]
    path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-wsrc-{inv_type_name}/traj-{traj}"
    if get_load_path(path_s + "/geon-info.txt") is not None:
        assert get_load_path(path_sp + "/checkpoint.txt") is not None
        q.displayln_info(f"{fname}: {job_tag}/{traj} already done.")
        return None
    if not check_prop_wsrc_avail_from_cps(job_tag, traj, inv_type=inv_type):
        q.displayln_info(f"{fname}: {job_tag}/{traj} {inv_type_name} CPS data not available.")
        return None
    q.check_stop()
    q.check_time_limit()
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return None
    with q.TimerFork(max_call_times_for_always_show_info=0, verbose=1, show_display=True):
        fsel = get_fsel()
        psel = get_psel()
        sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
        qar_sp = q.open_qar_info(get_save_path(path_sp + ".qar"), "a")
        total_site = q.Coordinate(get_param(job_tag, "total_site"))
        geo = q.Geometry(total_site)
        for inv_acc in [ 0, 1, 2, ]:
            for tslice in range(total_site[3]):
                fn = get_fn_prop_wsrc_avail_from_cps(
                    job_tag, traj,
                    tslice=tslice, inv_type=inv_type, inv_acc=inv_acc,
                    )
                if fn is None:
                    continue
                @q.timer_verbose
                def load_prop():
                    q.check_stop()
                    q.check_time_limit()
                    prop = q.Prop(geo)
                    qc.load_cps_prop_double(prop, fn)
                    return prop
                save_prop_wsrc_sparse(
                    job_tag, traj,
                    load_prop=load_prop,
                    tslice=tslice, inv_type=inv_type, inv_acc=inv_acc,
                    sfw=sfw, qar_sp=qar_sp,
                    psel=psel, fsel=fsel,
                    )
        sfw.close()
        qar_sp.write("checkpoint.txt", "", "", skip_if_exist=True)
        qar_sp.flush()
        qar_sp.close()
        q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
        is_doing_something = True
        q.release_lock()
    return is_doing_something

@q.timer_verbose
def run_job(job_tag, traj):
    fname = q.get_fname()
    if not check_prop_wsrc_avail_from_cps(job_tag, traj):
        q.displayln_info(f"{fname}: {job_tag}/{traj} CPS data not available.")
        return None
    is_doing_something = None
    with q.TimerFork(max_call_times_for_always_show_info=0, verbose=1, show_display=True):
        b = run_gauge_transform_from_cps(job_tag, traj)
        if b:
            is_doing_something = True
        with q.TimerFork(max_call_times_for_always_show_info=0, verbose=1, show_display=True):
            get_fsel_prob = run_fsel_prob_uniform(job_tag, traj)
            get_psel_prob = run_psel_prob_uniform(job_tag, traj)
            get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
            get_psel = run_psel_from_psel_prob(get_psel_prob)
        for inv_type in [ 0, 1, 2, ]:
            b = run_prop_wsrc_sparse_from_cps(
                job_tag, traj,
                inv_type=inv_type,
                get_fsel=get_fsel,
                get_psel=get_psel,
            )
            if b:
                is_doing_something = True
    return is_doing_something