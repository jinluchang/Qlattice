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
            assert False
    get_gf = ru.load_config_lazy(job_tag, path_gf)
    return get_gf

@q.timer_verbose
def run_gt(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    path_gt = get_load_path(f"gauge-transform/{job_tag}/traj={traj}.field")
    if path_gt is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-gauge_fix_coulomb"):
            gf = get_gf()
            gt = qg.gauge_fix_coulomb(gf)
            gt.save_double(get_save_path(f"gauge-transform/{job_tag}/traj={traj}.field"))
            q.release_lock()
            return lambda : gt
        else:
            return None
    else:
        @q.timer_verbose
        def load_gt():
            gt = q.GaugeTransform()
            gt.load_double(path_gt)
            # ADJUST ME
            # qg.check_gauge_fix_coulomb(get_gf(), gt)
            #
            return gt
        get_gt = q.lazy_call(load_gt)
    return get_gt

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
    path_psel = get_load_path(f"point-selection/{job_tag}/traj={traj}.txt")
    if path_psel is None:
        psel = mk_rand_psel(job_tag, traj)
        psel.save(get_save_path(f"point-selection/{job_tag}/traj={traj}.txt"))
        return lambda : psel
    else:
        @q.timer_verbose
        def load_psel():
            psel = q.PointSelection()
            psel.load(path_psel)
            return psel
        return q.lazy_call(load_psel)

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

@q.timer
def load_wall_src_info(path):
    # wi is a list of [ idx tslice inv_type inv_acc ]
    dt = q.qload_datatable_sync_node(path, True)
    t = [ list(map(int, l)) for l in dt ]
    wi = [ [ l[0], l[1], l[2], l[3], ] for l in t ]
    return wi

@q.timer_verbose
def run_wi(job_tag, traj):
    path_light = get_load_path(f"wall-src-info-light/{job_tag}/traj={traj}.txt")
    if path_light is None:
        wi_light = mk_rand_wall_src_info(job_tag, traj, inv_type = 0)
        save_wall_src_info(wi_light, get_save_path(f"wall-src-info-light/{job_tag}/traj={traj}.txt"));
    path_strange = get_load_path(f"wall-src-info-strange/{job_tag}/traj={traj}.txt")
    if path_strange is None:
        wi_strange = mk_rand_wall_src_info(job_tag, traj, inv_type = 1)
        save_wall_src_info(wi_strange, get_save_path(f"wall-src-info-strange/{job_tag}/traj={traj}.txt"));
    @q.timer_verbose
    def load():
        wi_light = load_wall_src_info(get_load_path(f"wall-src-info-light/{job_tag}/traj={traj}.txt"))
        wi_strange = load_wall_src_info(get_load_path(f"wall-src-info-strange/{job_tag}/traj={traj}.txt"))
        return wi_light + wi_strange
    return q.lazy_call(load)

@q.timer_verbose
def run_pi(job_tag, traj, get_psel):
    path = get_load_path(f"point-src-info/{job_tag}/traj={traj}.txt")
    if path is None:
        pi = mk_rand_point_src_info(job_tag, traj, get_psel())
        save_point_src_info(pi, get_save_path(f"point-src-info/{job_tag}/traj={traj}.txt"));
    @q.timer_verbose
    def load():
        pi = load_point_src_info(get_load_path(f"point-src-info/{job_tag}/traj={traj}.txt"))
        return pi
    return q.lazy_call(load)
