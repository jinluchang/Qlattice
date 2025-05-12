import numpy as np
import qlat as q
import qlat_gpt as qg
import gpt as g

from qlat_scripts.v1 import *

def set_param_field_selection_rate(job_tag, f_rate, n_points):
    """
    Use f_rate and n_points with Christoph & Mattia's strategy.
    #
    e.g.
    #
    set_param_field_selection_rate("96I", 1 / 64, 4096)
    """
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_volume = total_site.volume()
    n_points_fsel = round(total_volume * f_rate)
    fsel_prob = 1 - (1 - 1 / total_volume)**n_points_fsel
    psel_prob = 1 - (1 - 1 / total_volume)**n_points
    set_param(job_tag, "field_selection_fsel_rate")(fsel_prob)
    set_param(job_tag, "field_selection_psel_rate")(psel_prob)

def get_num_points_fsel_sampling(job_tag):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_volume = total_site.volume()
    fsel_rate = get_param(job_tag, "field_selection_fsel_rate")
    num_points_fsel_sampling = round(np.log(1 - fsel_rate) / np.log(1 - 1 / total_volume))
    return num_points_fsel_sampling

def get_num_points_psel_sampling(job_tag):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_volume = total_site.volume()
    psel_rate = get_param(job_tag, "field_selection_psel_rate")
    num_points_psel_sampling = round(np.log(1 - psel_rate) / np.log(1 - 1 / total_volume))
    return num_points_psel_sampling

@q.cache_call()
@q.timer_verbose
def get_all_positions(job_tag, traj):
    """
    Based on the email:
    #
    From: Christoph Lehner <christoph@lhnr.de>
    Date: Thu, 5 Dec 2024 09:50:26 +0100
    Subject: Re: 96I propagators
    To: chulwoo@bnl.gov
    Cc: "Jin, Luchang" <luchang.jin@uconn.edu>, Luchang Jin <ljin.luchang@gmail.com>
    """
    assert isinstance(traj, int)
    conf = traj
    rng = g.random(f"{job_tag}-positions-{conf}")
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    total_volume = total_site.volume()
    num_points_fsel_sampling = get_num_points_fsel_sampling(job_tag)
    n_all_points = num_points_fsel_sampling
    # n_all_points = 2654208 # 1/64 of all points, 96I
    L = total_site.to_list()
    all_positions = np.array([
        [rng.uniform_int(min=0, max=L[i] - 1) for i in range(4)] for j in range(n_all_points)
    ], dtype=np.int32)
    return all_positions

@q.cache_call()
@q.timer_verbose
def get_psrc_positions(job_tag, traj):
    """
    Based on the email:
    #
    From: Christoph Lehner <christoph@lhnr.de>
    Date: Fri, 6 Dec 2024 17:30:43 +0100
    Subject: Re: 96I propagators
    To: Luchang Jin <ljin.luchang@gmail.com>
    Cc: chulwoo@bnl.gov
    """
    num_points_psel_sampling = get_num_points_psel_sampling(job_tag)
    all_positions = get_all_positions(job_tag, traj)
    psrc_positions = all_positions[:num_points_psel_sampling].copy()
    return psrc_positions

@q.timer
def run_fsel_prob_uniform(job_tag, traj):
    """
    return get_fsel_prob
        fsel_prob = get_fsel_prob()
        fsel = fsel_prob.fsel
    #
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = lambda : get_fsel_prob().fsel
    #
    `fsel` for 96I from Christoph and Mattia.
    """
    fname = q.get_fname()
    fn_fsel = f"{job_tag}/field-selection/traj-{traj}.field"
    fn_fsel_prob = f"{job_tag}/field-selection-weight/traj-{traj}/fsel-prob.sfield"
    @q.lazy_call
    @q.timer_verbose
    def get_fsel_prob():
        fsel = q.FieldSelection()
        total_size = fsel.load(get_load_path(fn_fsel))
        assert total_size > 0
        fsel_prob = q.SelectedFieldRealD(fsel, 1)
        total_size = fsel_prob.load_double(get_load_path(fn_fsel_prob))
        assert total_size > 0
        return fsel_prob
    ret = get_fsel_prob
    if get_load_path(fn_fsel) is not None:
        if get_load_path(fn_fsel_prob) is not None:
            return ret
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    fsel_rate = get_param(job_tag, "field_selection_fsel_rate")
    q.displayln_info(-1, fname, f"fsel_rate = {fsel_rate}")
    assert fsel_rate is not None
    assert get_load_path(fn_fsel) is None
    assert get_load_path(fn_fsel_prob) is None
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    xg_arr = get_all_positions(job_tag, traj)
    psel = q.PointsSelection(total_site, xg_arr)
    fsel = q.FieldSelection(psel)
    fsel.save(get_save_path(fn_fsel))
    fsel_prob = q.SelectedFieldRealD(fsel, 1)
    fsel_prob[:] = fsel_rate
    fsel_prob.save_double(get_save_path(fn_fsel_prob))
    q.release_lock()
    return ret

@q.timer
def run_psel_prob_uniform(job_tag, traj):
    """
    return get_psel_prob
        psel_prob = get_psel_prob()
        psel = psel_prob.psel
    #
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel = lambda : get_psel_prob().psel
    """
    fname = q.get_fname()
    fn_psel = f"{job_tag}/point-selection/traj-{traj}.txt"
    fn_psel_prob = f"{job_tag}/field-selection-weight/traj-{traj}/psel-prob.lat"
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    @q.lazy_call
    @q.timer_verbose
    def get_psel_prob():
        psel = q.PointsSelection()
        psel.load(get_load_path(fn_psel), geo)
        psel_prob = q.SelectedPointsRealD(psel, 1)
        psel_prob.load(get_load_path(fn_psel_prob))
        return psel_prob
    ret = get_psel_prob
    if get_load_path(fn_psel) is not None:
        if get_load_path(fn_psel_prob) is not None:
            return ret
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    psel_rate = get_param(job_tag, "field_selection_psel_rate")
    q.displayln_info(-1, fname, f"psel_rate = {psel_rate}")
    assert psel_rate is not None
    assert get_load_path(fn_psel) is None
    assert get_load_path(fn_psel_prob) is None
    xg_arr = get_all_positions(job_tag, traj)
    psel = q.PointsSelection(total_site, xg_arr)
    fsel = q.FieldSelection(psel)
    psel = fsel.to_psel()
    psel.save(get_save_path(fn_psel))
    psel_prob = q.SelectedPointsRealD(psel, 1)
    psel_prob[:] = psel_rate
    psel_prob.save(get_save_path(fn_psel_prob))
    q.release_lock()
    return ret