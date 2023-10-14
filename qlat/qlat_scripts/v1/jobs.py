import qlat as q
from . import rbc_ukqcd_params as rup
from .rbc_ukqcd_params import set_param, get_param

import numpy as np

import pprint
import functools
import os

save_path_default = "results"

load_path_list = [ "results", ]

def get_save_path(fn):
    return os.path.join(save_path_default, fn)

def get_load_path(*fns):
    def get(fn):
        if fn is None:
            return None
        elif isinstance(fn, (tuple, list)):
            for f in fn:
                p = get(f)
                if p is not None:
                    return p
        else:
            for path in load_path_list:
                p = os.path.join(path, fn)
                if q.does_file_exist_qar_sync_node(p):
                    return p
        return None
    return get(fns)

# ----------

@q.timer_verbose
def check_job(job_tag, traj, fns_produce, fns_need):
    """
    return False if config is finished or unavailable
    """
    is_job_done = True
    for fn in fns_produce:
        if get_load_path(fn) is None:
            q.displayln_info(f"check_job: {job_tag} {traj} to do as '{fn}' does not exist.")
            is_job_done = False
            break
    if is_job_done:
        return False
    #
    is_job_avail = True
    for fn in fns_need:
        if get_load_path(fn) is None:
            q.displayln_info(f"check_job: {job_tag} {traj} unavailable as '{fn}' does not exist.")
            is_job_avail = False
            break
    if not is_job_avail:
        return False
    #
    q.check_stop()
    q.check_time_limit()
    #
    assert not is_job_done and is_job_avail
    #
    return True

# ----------

@q.timer_verbose
def run_gf(job_tag, traj):
    path_gf = get_load_path(
            f"{job_tag}/configs/ckpoint_lat.{traj}",
            f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",
            )
    if path_gf is None:
        if job_tag[:5] == "test-":
            if not q.obtain_lock(f"locks/{job_tag}-{traj}-gauge_field"):
                return None
            total_site = rup.get_total_site(job_tag)
            gf = rup.mk_sample_gauge_field_v3(job_tag, f"{traj}")
            path_gf = get_save_path(f"{job_tag}/configs/ckpoint_lat.{traj}")
            gf.save(path_gf)
            q.release_lock()
        else:
            @q.timer_verbose
            def load_gf():
                assert False
            return load_gf
    get_gf = rup.load_config_lazy(job_tag, path_gf)
    return get_gf

@q.timer_verbose
def run_gt(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    tfn = f"{job_tag}/gauge-transform/traj-{traj}.field"
    path_gt = get_load_path(tfn)
    if path_gt is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-gauge_fix_coulomb"):
            gf = get_gf()
            import qlat_gpt as qg
            gt = qg.gauge_fix_coulomb(gf)
            gt.save_cps(get_save_path(f"{job_tag}/gauge-transform/traj-{traj}.gfix"))
            gt.save_double(get_save_path(tfn))
            q.release_lock()
            return lambda : gt
        else:
            return None
    @q.timer_verbose
    def load_gt():
        path_gt = get_load_path(tfn)
        assert path_gt is not None
        gt = q.GaugeTransform()
        gt.load_double(path_gt)
        # ADJUST ME
        # import qlat_gpt as qg
        # qg.check_gauge_fix_coulomb(get_gf(), gt)
        #
        return gt
    get_gt = q.lazy_call(load_gt)
    return get_gt

# ----------

@q.timer
def mk_rand_wall_src_info_n_exact(job_tag, traj, inv_type):
    params = rup.dict_params[job_tag]
    n_exact = params["n_exact_wsrc"]
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_wall_src_info")
    inv_acc_s = 1
    inv_acc_e = 2
    total_site = rup.get_total_site(job_tag)
    t_size = total_site[3]
    wi_s = [ [ t, inv_type, inv_acc_s, ] for t in range(t_size) ]
    mask = [ False, ] * t_size
    for i in range(n_exact):
        t_e = rs.rand_gen() % t_size
        mask[t_e] = True
    wi_e = []
    for t in range(t_size):
        if mask[t]:
            wi_e.append([ t, inv_type, inv_acc_e, ])
    wi = wi_e + wi_s
    for i in range(len(wi)):
        wi[i] = [ i, ] + wi[i]
    return wi

@q.timer
def mk_rand_wall_src_info_prob(job_tag, traj, inv_type):
    params = rup.dict_params[job_tag]
    prob = params["prob_exact_wsrc"]
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_wall_src_info_prob")
    inv_acc_s = 1
    inv_acc_e = 2
    total_site = rup.get_total_site(job_tag)
    t_size = total_site[3]
    wi_s = [ [ t, inv_type, inv_acc_s, ] for t in range(t_size) ]
    wi_e = []
    for t in range(t_size):
        if rs.u_rand_gen() < prob:
            wi_e.append([ t, inv_type, inv_acc_e, ])
    wi = wi_e + wi_s
    for i in range(len(wi)):
        wi[i] = [ i, ] + wi[i]
    return wi

def get_prob_exact_wsrc(job_tag):
    params = rup.dict_params[job_tag]
    if "prob_exact_wsrc" in params:
        return params["prob_exact_wsrc"]
    n_exact = params["n_exact_wsrc"]
    total_site = rup.get_total_site(job_tag)
    return 1 - (1 - 1 / total_site[3])**n_exact

@q.timer
def mk_rand_wall_src_info(job_tag, traj, inv_type):
    """
    wi is a list of [ idx tslice inv_type inv_acc ]
    """
    params = rup.dict_params[job_tag]
    if "prob_exact_wsrc" not in params:
        return mk_rand_wall_src_info_n_exact(job_tag, traj, inv_type)
    return mk_rand_wall_src_info_prob(job_tag, traj, inv_type)

@q.timer
def save_wall_src_info(wi, path):
    """
    wi is a list of  [ idx tslice inv_type inv_acc ]
    """
    if 0 != q.get_id_node():
        return None
    lines = [ " ".join([ f"{v:5d}" for v in l ]) for l in wi ]
    content = "\n".join(lines + [ "", ])
    q.qtouch(path, content)

@q.timer
def load_wall_src_info(path):
    assert path is not None
    """
    wi is a list of [ idx tslice inv_type inv_acc ]
    """
    dt = q.qload_datatable_sync_node(path, True)
    t = [ list(map(int, l)) for l in dt ]
    wi = [ [ l[0], l[1], l[2], l[3], ] for l in t ]
    return wi

@q.timer_verbose
def run_wi(job_tag, traj):
    tfn_l = f"{job_tag}/wall-src-info-light/traj-{traj}.txt"
    tfn_s = f"{job_tag}/wall-src-info-strange/traj-{traj}.txt"
    path_light = get_load_path(tfn_l)
    if path_light is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-wi"):
            wi_light = mk_rand_wall_src_info(job_tag, traj, inv_type = 0)
            save_wall_src_info(wi_light, get_save_path(tfn_l));
            q.release_lock()
        else:
            return None
    path_strange = get_load_path(tfn_s)
    if path_strange is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-wi"):
            wi_strange = mk_rand_wall_src_info(job_tag, traj, inv_type = 1)
            save_wall_src_info(wi_strange, get_save_path(tfn_s));
            q.release_lock()
        else:
            return None
    @q.timer_verbose
    def load():
        wi_light = load_wall_src_info(get_load_path(tfn_l))
        wi_strange = load_wall_src_info(get_load_path(tfn_s))
        return wi_light + wi_strange
    return q.lazy_call(load)

# ----------

def get_n_points_psel(job_tag):
    assert job_tag in rup.dict_params
    assert "n_points_psel" in rup.dict_params[job_tag]
    return get_param(job_tag, "n_points_psel")

@q.timer
def mk_rand_psel(job_tag, traj):
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_psel")
    total_site = rup.get_total_site(job_tag)
    n_points = get_n_points_psel(job_tag)
    psel = q.PointsSelection()
    psel.set_rand(rs, total_site, n_points)
    return psel

@q.timer_verbose
def run_psel(job_tag, traj):
    tfn = f"{job_tag}/point-selection/traj-{traj}.txt"
    path_psel = get_load_path(tfn)
    if path_psel is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-psel"):
            psel = mk_rand_psel(job_tag, traj)
            psel.save(get_save_path(tfn))
            q.release_lock()
        else:
            return None
    #
    @q.timer_verbose
    def load_psel():
        path_psel = get_load_path(tfn)
        assert path_psel is not None
        total_site = rup.get_total_site(job_tag)
        psel = q.PointsSelection()
        psel.load(path_psel, q.Geometry(total_site))
        assert psel.n_points() == get_n_points_psel(job_tag)
        return psel
    return q.lazy_call(load_psel)

# ----------

def get_n_points_pi(job_tag, traj, inv_type, inv_acc):
    assert job_tag in rup.dict_params
    assert "n_points" in rup.dict_params[job_tag]
    return rup.dict_params[job_tag]["n_points"][inv_type][inv_acc]

@q.timer
def mk_rand_point_src_info(job_tag, traj, psel):
    # pi is a list of [ idx xg inv_type inv_acc ]
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_point_src_info")
    xg_list = psel.to_list()
    assert len(xg_list) == get_n_points_pi(job_tag, traj, 0, 0)
    g_pi = [ [] for _ in xg_list ]
    for inv_type in [ 0, 1, ]:
        for inv_acc in [ 0, 1, 2, ]:
            for i in range(get_n_points_pi(job_tag, traj, inv_type, inv_acc)):
                g_pi[i].append([ xg_list[i], inv_type, inv_acc ])
    pi = []
    for g in g_pi:
        pi += g
    for i in range(len(pi)):
        pi[i] = [ i, ] + pi[i]
    return pi

@q.timer
def save_point_src_info(pi, path):
    # pi is a list of [ idx xg inv_type inv_acc ]
    if 0 != q.get_id_node():
        return None
    def mk_line(l):
        [ idx, xg, inv_type, inv_acc ] = l
        return f"{idx:5d}    {xg[0]:3d} {xg[1]:3d} {xg[2]:3d} {xg[3]:3d}    {inv_type:3d} {inv_acc:3d}"
    lines = list(map(mk_line, pi))
    content = "\n".join([ f"{len(lines)}" ] + lines + [ "" ])
    q.qtouch(path, content)

@q.timer
def load_point_src_info(path):
    # pi is a list of [ idx xg inv_type inv_acc ]
    dt = q.qload_datatable_sync_node(path, True)
    t = [ list(map(int, l)) for l in dt ][1:]
    pi = [ [ l[0], l[1:5], l[5], l[6], ] for l in t ]
    return pi

@q.timer_verbose
def run_pi(job_tag, traj, get_psel):
    tfn = f"{job_tag}/point-src-info/traj-{traj}.txt"
    path = get_load_path(tfn)
    if path is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-pi"):
            pi = mk_rand_point_src_info(job_tag, traj, get_psel())
            save_point_src_info(pi, get_save_path(tfn));
            q.release_lock()
        else:
            return None
    @q.timer_verbose
    def load():
        path = get_load_path(tfn)
        assert path is not None
        pi = load_point_src_info(path)
        return pi
    return q.lazy_call(load)

# ----------

@q.timer
def load_point_distribution(job_tag):
    """
    return point_distribution
    where
    point_distribution[xg_rel] = probability of the relative coordinate of a point relative to a selected point equals to ``xg_rel``.
    xg_rel = (x, y, z, t,)
    x >= y >= z >= 0 and t >= 0
    n_points = get_n_points_psel(job_tag)
    """
    n_points = get_n_points_psel(job_tag)
    tfn = f"{job_tag}/point-distribution/point-distribution.txt"
    path = get_load_path(tfn)
    if path is None:
        return None
    dt = q.qload_datatable_sync_node(path, True)
    point_distribution = dict()
    for l in dt:
        x, y, z, t, prob = l
        x = int(x)
        y = int(y)
        z = int(z)
        t = int(t)
        point_distribution[(x, y, z, t,)] = prob * (n_points - 1) / n_points
    point_distribution[(0, 0, 0, 0,)] = 1.0 / n_points
    return point_distribution

def classify_rel_coordinate(xg_rel_arrary, total_site_array):
    """
    xg_rel_arrary = np.array(xg_rel)
    total_site_array = np.array(total_site)
    """
    total_site_half = total_site_array // 2
    xg_rel_arrary = xg_rel_arrary % total_site_array
    xg_rel_abs = total_site_half - abs(xg_rel_arrary - total_site_half)
    x, y, z, t = xg_rel_abs
    x, y, z = sorted([x, y, z])
    return (x, y, z, t,)

def get_point_xrel_prob(xg_rel_arrary, total_site_array, point_distribution, n_points):
    """
    xg_rel_arrary = np.array(xg_rel)
    total_site_array = np.array(total_site)
    point_distribution = load_point_distribution(job_tag)
    n_points = get_n_points_psel(job_tag)
    """
    if point_distribution is None:
        if np.all(xg_rel_arrary == 0):
            return 1.0 / n_points
        else:
            total_volume = np.prod(total_site_array)
            return (n_points - 1) / n_points / (total_volume - 1)
    xg_rel = classify_rel_coordinate(xg_rel_arrary, total_site_array)
    prob = point_distribution[xg_rel]
    return prob

# ----------

@q.timer
def mk_rand_fsel(job_tag, traj, n_per_tslice):
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_fsel")
    total_site = rup.get_total_site(job_tag)
    fsel = q.FieldSelection()
    fsel.set_rand(rs, total_site, n_per_tslice)
    return fsel

@q.timer
def mk_fselc(fsel, psel):
    fselc = fsel.copy()
    fselc.add_psel(psel)
    return fselc

@q.timer_verbose
def run_fsel(job_tag, traj, get_psel):
    if get_psel is None:
        return None
    tfn = f"{job_tag}/field-selection/traj-{traj}.field"
    path_fsel = get_load_path(tfn)
    total_site = rup.get_total_site(job_tag)
    n_per_tslice = total_site[0] * total_site[1] * total_site[2] // 16
    if path_fsel is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-fsel"):
            fsel = mk_rand_fsel(job_tag, traj, n_per_tslice)
            fsel.save(get_save_path(tfn))
            fselc = mk_fselc(fsel, get_psel())
            q.release_lock()
            return lambda : ( fsel, fselc, )
        else:
            return None
    @q.timer_verbose
    def load_fsel():
        path_fsel = get_load_path(tfn)
        assert path_fsel is not None
        fsel = q.FieldSelection()
        total_size = fsel.load(path_fsel, n_per_tslice)
        assert total_size > 0
        fselc = mk_fselc(fsel, get_psel())
        return fsel, fselc
    return q.lazy_call(load_fsel)

# ----------

@q.timer
def mk_rand_fsel_smear(job_tag, traj, n_per_tslice_smear):
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_fsel_smear")
    total_site = rup.get_total_site(job_tag)
    fsel = q.FieldSelection()
    fsel.set_rand(rs, total_site, n_per_tslice_smear)
    return fsel

@q.timer_verbose
def run_psel_smear(job_tag, traj):
    # return lambda : psel_smear
    # psel_smear should randomly select same number of point on each tslice
    tfn = f"{job_tag}/point-selection-smear/traj-{traj}.txt"
    path_psel = get_load_path(tfn)
    total_site = rup.get_total_site(job_tag)
    if path_psel is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-psel-smear"):
            n_per_tslice_smear = rup.dict_params[job_tag]["n_per_tslice_smear"]
            fsel = mk_rand_fsel_smear(job_tag, traj, n_per_tslice_smear)
            psel = fsel.to_psel()
            psel.save(get_save_path(tfn))
            q.release_lock()
        else:
            return None
    #
    @q.timer_verbose
    def load_psel():
        path_psel = get_load_path(tfn)
        assert path_psel is not None
        total_site = rup.get_total_site(job_tag)
        psel = q.PointsSelection()
        psel.load(path_psel, q.Geometry(total_site))
        return psel
    return q.lazy_call(load_psel)

# ----------

@q.timer
def run_gf_ape(job_tag, get_gf):
    if get_gf is None:
        return None
    coef = rup.dict_params[job_tag]["gf_ape_smear_coef"]
    step = rup.dict_params[job_tag]["gf_ape_smear_step"]
    #
    @q.timer_verbose
    def run():
        gf = get_gf()
        gf_ape = gf.copy()
        q.gf_spatial_ape_smear(gf_ape, coef, step)
        gf_ape = q.mk_left_expanded_gauge_field(gf_ape)
        return gf_ape
    return q.lazy_call(run)

# ----------

@q.timer
def compute_eig(gf, job_tag, inv_type = 0, inv_acc = 0, *, path = None):
    """
    return a function ``get_eig''
    ``get_eig()'' return the ``eig''
    """
    from . import rbc_ukqcd as ru
    load_eig = ru.load_eig_lazy(get_load_path(path), job_tag)
    if load_eig is not None:
        return load_eig
    # evec, evals = ru.mk_eig(gf, job_tag, inv_type, inv_acc)
    basis, cevec, smoothed_evals = ru.mk_ceig(gf, job_tag, inv_type, inv_acc)
    eig = [ basis, cevec, smoothed_evals, ]
    ru.save_ceig(get_save_path(path + ".partial"), eig, job_tag, inv_type, inv_acc);
    q.qrename_info(get_save_path(path + ".partial"), get_save_path(path))
    test_eig(gf, eig, job_tag, inv_type)
    def get_eig():
        return eig
    return get_eig

@q.timer
def test_eig(gf, eig, job_tag, inv_type):
    from . import rbc_ukqcd as ru
    geo = gf.geo()
    src = q.FermionField4d(geo)
    src.set_rand(q.RngState("test_eig:src.set_rand"))
    q.displayln_info(f"src norm {src.qnorm():.10E}")
    sol_ref = ru.get_inv(gf, job_tag, inv_type, inv_acc = 2, eig = eig, eps = 1e-10, mpi_split = False, qtimer = False) * src
    q.displayln_info(f"sol_ref norm {sol_ref.qnorm():.10E} with eig")
    for inv_acc in [ 0, 1, 2, ]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig, mpi_split = False, qtimer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} with eig")
        if inv_acc in [ 0, 1, ]:
            sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, mpi_split = False, qtimer = False) * src
            sol -= sol_ref
            q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} without eig")

@q.timer_verbose
def run_eig(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    from . import rbc_ukqcd as ru
    get_eig = ru.load_eig_lazy(get_load_path(f"{job_tag}/eig/traj-{traj}"), job_tag)
    if get_eig is None and get_gf is not None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-run-eig"):
            get_eig = compute_eig(get_gf(), job_tag, inv_type = 0, path = f"{job_tag}/eig/traj-{traj}")
            q.release_lock()
            return get_eig
        else:
            return None
    else:
        return get_eig

@q.timer_verbose
def run_eig_strange(job_tag, traj, get_gf):
    # if failed, return None
    # if no parameter, return lambda : None
    if None in [ get_gf, ]:
        return None
    if 1 not in rup.dict_params[job_tag]["clanc_params"]:
        fn = f"{job_tag}/eig-strange/traj-{traj}/no-eig-parameters.txt"
        if get_load_path(fn) is None:
            q.qtouch_info(get_save_path(fn))
        return lambda : None
    from . import rbc_ukqcd as ru
    get_eig = ru.load_eig_lazy(get_load_path(f"{job_tag}/eig-strange/traj-{traj}"), job_tag)
    if get_eig is None and get_gf is not None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-run-eig-strange"):
            get_eig = compute_eig(get_gf(), job_tag, inv_type = 1, path = f"{job_tag}/eig-strange/traj-{traj}")
            q.release_lock()
            return get_eig
        else:
            return None
    else:
        return get_eig

# ----------

@functools.lru_cache(maxsize=None)
def get_r_list(job_tag):
    total_site = q.Coordinate(rup.dict_params[job_tag]["total_site"])
    r_limit = q.get_r_limit(total_site)
    r_list = q.mk_r_list(r_limit, r_all_limit=28.0, r_scaling_factor=5.0)
    # r_list = q.mk_r_list(r_limit, r_all_limit=0.0, r_scaling_factor=5.0) # old choice
    return r_list

@functools.lru_cache(maxsize=None)
def get_r_sq_interp_idx_coef_list(job_tag):
    """
    Return [ (r_idx_low, r_idx_high, coef_low, coef_high,), ... ] indexed by r_sq
    """
    r_list = get_r_list(job_tag)
    return q.mk_r_sq_interp_idx_coef_list(r_list)

@q.timer_verbose
def run_r_list(job_tag):
    fn = f"{job_tag}/r_list/r_list.lat"
    r_list = get_r_list(job_tag)
    ld = q.mk_lat_data([
        [ "r_idx", len(r_list), ],
        ])
    ld.from_numpy(np.array(r_list))
    if get_load_path(fn) is not None:
        ld_load = q.load_lat_data(get_load_path(fn))
        assert ld.is_match(ld_load)
        ld_diff = ld - ld_load
        assert ld_diff.qnorm() < 1e-20
        return
    ld.save(get_save_path(fn))

# ----------
