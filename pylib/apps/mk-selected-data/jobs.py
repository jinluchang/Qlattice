import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import os

save_path_default = "results"

load_path_list = [ "results", ]

def get_save_path(fn):
    return os.path.join(save_path_default, fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = load_path_list
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

# ----------

@q.timer_verbose
def run_gf(job_tag, traj):
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    if path_gf is None:
        if job_tag[:5] == "test-":
            gf = ru.mk_sample_gauge_field(job_tag, f"{traj}")
            path_gf = get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
            # gf.save(path_gf)
            qg.save_gauge_field(gf, path_gf)
        else:
            @q.timer_verbose
            def load_gf():
                assert False
            return load_gf
    get_gf = ru.load_config_lazy(job_tag, path_gf)
    return get_gf

@q.timer_verbose
def run_gt(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    tfn = f"gauge-transform/{job_tag}/traj={traj}.field"
    path_gt = get_load_path(tfn)
    if path_gt is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-gauge_fix_coulomb"):
            gf = get_gf()
            gt = qg.gauge_fix_coulomb(gf)
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
        # qg.check_gauge_fix_coulomb(get_gf(), gt)
        #
        return gt
    get_gt = q.lazy_call(load_gt)
    return get_gt

# ----------

@q.timer
def mk_rand_wall_src_info(job_tag, traj, inv_type):
    # wi is a list of [ idx tslice inv_type inv_acc ]
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_wall_src_info")
    inv_acc_s = 1
    inv_acc_e = 2
    total_site = ru.get_total_site(job_tag)
    t_size = total_site[3]
    wi_s = [ [ t, inv_type, inv_acc_s ] for t in range(t_size) ]
    n_exact = 2
    mask = [ False, ] * t_size
    for i in range(n_exact):
        t_e = rs.rand_gen() % t_size
        mask[t_e] = True
    wi_e = []
    for t in range(t_size):
        if mask[t]:
            wi_e.append([ t, inv_type, inv_acc_e ])
    wi = wi_e + wi_s
    for i in range(len(wi)):
        wi[i] = [ i, ] + wi[i]
    return wi

@q.timer
def save_wall_src_info(wi, path):
    # wi is a list of  [ idx tslice inv_type inv_acc ]
    if 0 != q.get_id_node():
        return None
    lines = [ " ".join([ f"{v:5d}" for v in l ]) for l in wi ]
    content = "\n".join(lines + [ "", ])
    q.qtouch(path, content)

@q.timer
def load_wall_src_info(path):
    assert path is not None
    # wi is a list of [ idx tslice inv_type inv_acc ]
    dt = q.qload_datatable_sync_node(path, True)
    t = [ list(map(int, l)) for l in dt ]
    wi = [ [ l[0], l[1], l[2], l[3], ] for l in t ]
    return wi

@q.timer_verbose
def run_wi(job_tag, traj):
    tfn_l = f"wall-src-info-light/{job_tag}/traj={traj}.txt"
    tfn_s = f"wall-src-info-strange/{job_tag}/traj={traj}.txt"
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

def get_n_points(job_tag, traj, inv_type, inv_acc):
    assert job_tag in rup.dict_params
    assert "n_points" in rup.dict_params[job_tag]
    return rup.dict_params[job_tag]["n_points"][inv_type][inv_acc]

@q.timer
def mk_rand_psel(job_tag, traj):
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_psel")
    total_site = ru.get_total_site(job_tag)
    n_points = get_n_points(job_tag, traj, 0, 0)
    psel = q.PointSelection()
    psel.set_rand(rs, total_site, n_points)
    return psel

@q.timer_verbose
def run_psel(job_tag, traj):
    tfn = f"point-selection/{job_tag}/traj={traj}.txt"
    path_psel = get_load_path(tfn)
    if path_psel is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-psel"):
            psel = mk_rand_psel(job_tag, traj)
            psel.save(get_save_path(tfn))
            q.release_lock()
            return lambda : psel
        else:
            return None
    @q.timer_verbose
    def load_psel():
        path_psel = get_load_path(tfn)
        assert path_psel is not None
        psel = q.PointSelection()
        psel.load(path_psel)
        return psel
    return q.lazy_call(load_psel)

# ----------

@q.timer
def mk_rand_point_src_info(job_tag, traj, psel):
    # pi is a list of [ idx xg inv_type inv_acc ]
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_point_src_info")
    xg_list = psel.to_list()
    assert len(xg_list) == get_n_points(job_tag, traj, 0, 0)
    g_pi = [ [] for _ in xg_list ]
    for inv_type in [ 0, 1, ]:
        for inv_acc in [ 0, 1, 2, ]:
            for i in range(get_n_points(job_tag, traj, inv_type, inv_acc)):
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
    tfn = f"point-src-info/{job_tag}/traj={traj}.txt"
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
def mk_rand_fsel(job_tag, traj, n_per_tslice):
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_fsel")
    total_site = ru.get_total_site(job_tag)
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
    tfn = f"field-selection/{job_tag}/traj={traj}.field"
    path_fsel = get_load_path(tfn)
    total_site = ru.get_total_site(job_tag)
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
def compute_eig(gf, job_tag, inv_type = 0, inv_acc = 0, *, path = None):
    # return a function ``get_eig''
    # ``get_eig()'' return the ``eig''
    load_eig = ru.load_eig_lazy(get_load_path(path), job_tag)
    if load_eig is not None:
        return load_eig
    # evec, evals = ru.mk_eig(gf, job_tag, inv_type, inv_acc)
    basis, cevec, smoothed_evals = ru.mk_ceig(gf, job_tag, inv_type, inv_acc)
    eig = [ basis, cevec, smoothed_evals ]
    ru.save_ceig(get_save_path(path + ".partial"), eig, job_tag, inv_type, inv_acc);
    q.qrename_info(get_save_path(path + ".partial"), get_save_path(path))
    test_eig(gf, eig, job_tag, inv_type)
    def get_eig():
        return eig
    return get_eig

@q.timer
def test_eig(gf, eig, job_tag, inv_type):
    geo = gf.geo()
    src = q.FermionField4d(geo)
    q.displayln_info(f"src norm {src.qnorm()}")
    src.set_rand(q.RngState("test_eig:{id(inv)}"))
    sol_ref = ru.get_inv(gf, job_tag, inv_type, inv_acc = 2, eig = eig, eps = 1e-10, mpi_split = False, timer = False) * src
    q.displayln_info(f"sol_ref norm {sol_ref.qnorm()} with eig")
    for inv_acc in [0, 1, 2]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig, mpi_split = False, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} with eig")
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, mpi_split = False, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} without eig")

@q.timer_verbose
def run_eig(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    get_eig = ru.load_eig_lazy(get_load_path(f"eig/{job_tag}/traj={traj}"), job_tag)
    if get_eig is None and get_gf is not None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-run-eig"):
            q.qmkdir_info(get_save_path(f"eig"))
            q.qmkdir_info(get_save_path(f"eig/{job_tag}"))
            get_eig = compute_eig(get_gf(), job_tag, inv_type = 0, path = f"eig/{job_tag}/traj={traj}")
            q.release_lock()
            return get_eig
        else:
            return None
    else:
        return get_eig

# ----------
