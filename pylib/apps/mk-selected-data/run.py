#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [
            "results",
            "../mk-gf-gt/results",
            "../mk-lanc/results",
            "/gpfs/alpine/lgt116/proj-shared/ljin",
            ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer_verbose
def check_job(job_tag, traj):
    # return True if config is finished or unavailable
    fns_produce = [
            get_load_path(f"point-selection/{job_tag}/traj={traj}.txt"),
            get_load_path(f"field-selection/{job_tag}/traj={traj}.field"),
            get_load_path(f"prop-wsrc-strange/{job_tag}/traj={traj}"),
            get_load_path(f"prop-wsrc-light/{job_tag}/traj={traj}"),
            get_load_path(f"prop-psrc-strange/{job_tag}/traj={traj}"),
            get_load_path(f"prop-psrc-light/{job_tag}/traj={traj}"),
            ]
    is_job_done = True
    for fn in fns_produce:
        if fn is None:
            q.displayln_info(f"check_job: {job_tag} {traj} to do as some file does not exist.")
            is_job_done = False
            break
    if is_job_done:
        return True
    #
    fns_need = [
            get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}"),
            get_load_path(f"gauge-transform/{job_tag}/traj={traj}.field"),
            get_load_path(f"eig/{job_tag}/traj={traj}"),
            ]
    for fn in fns_need:
        if fn is None:
            q.displayln_info(f"check_job: {job_tag} {traj} unavailable as {fn} does not exist.")
            return True
    #
    q.check_stop()
    q.check_time_limit()
    #
    q.qmkdir_info(f"locks")
    q.qmkdir_info(get_save_path(f""))
    #
    return False

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
    mask = [ False ] * t_size
    for i in range(n_exact):
        t_e = rs.rand_gen() % t_size
        mask[t_e] = True
    wi_e = []
    for t in range(t_size):
        if mask[t]:
            wi_e.append([ t, inv_type, inv_acc_e ])
    wi = wi_e + wi_s
    for i in range(len(wi)):
        wi[i] = [ i ] + wi[i]
    return wi

@q.timer
def save_wall_src_info(wi, path):
    # wi is a list of  [ idx tslice inv_type inv_acc ]
    if 0 != q.get_id_node():
        return None
    lines = [ " ".join([ f"{v:5d}" for v in l ]) for l in wi ]
    content = "\n".join(lines + [""])
    q.qtouch(get_save_path(path), content)

@q.timer
def mk_rand_point_src_info(job_tag, traj, psel):
    # pi is a list of [ idx xg inv_type inv_acc ]
    rs = q.RngState(f"seed {job_tag} {traj}").split("mk_rand_point_src_info")
    xg_list = psel.list()
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
        pi[i] = [ i ] + pi[i]
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
    q.qtouch(get_save_path(path), content)

@q.timer_verbose
def compute_prop(inv, src, *, tag, sfw, fn_sp, psel, fsel, fselc):
    sol = inv * src
    s_sol = sol.sparse(fselc)
    s_sol.save_float_from_double(sfw, tag)
    sp_sol = s_sol.sparse(psel)
    sp_sol.save(get_save_path(fn_sp))
    return sol

@q.timer
def compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc, *,
        idx, sfw, path_sp, psel, fsel, fselc, eig, finished_tags):
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_wsrc: idx={idx} tslice={tslice}", job_tag, inv_type, inv_acc)
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    src = q.mk_wall_src(geo, tslice)
    fn_sp = os.path.join(path_sp, f"{tag}.lat")
    prop = compute_prop(inv, src, tag = tag, sfw = sfw, fn_sp = fn_sp, psel = psel, fsel = fsel, fselc = fselc)
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    prop.glb_sum_tslice().save(get_save_path(fn_spw))

@q.timer_verbose
def compute_prop_wsrc_all(gf, gt, wi, job_tag, inv_type, *,
        path_s, path_sp, psel, fsel, fselc, eig):
    if q.does_file_exist_sync_node(get_save_path(path_s + ".acc.partial")):
        q.qrename_info(get_save_path(path_s + ".acc.partial"), get_save_path(path_s + ".acc"))
    finished_tags = q.properly_truncate_fields_sync_node(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 4, ])
    for inv_acc in [ 2, 1 ]:
        for p in wi:
            idx, tslice, inv_type_p, inv_acc_p = p
            if inv_type_p == inv_type and inv_acc_p == inv_acc:
                compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc,
                        idx = idx, sfw = sfw, path_sp = path_sp,
                        psel = psel, fsel = fsel, fselc = fselc, eig = eig,
                        finished_tags = finished_tags)
        q.displayln_info(pprint.pformat(q.list_cache()))
        q.clean_cache(q.cache_inv)
        q.displayln_info(pprint.pformat(q.list_cache()))
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

@q.timer
def compute_prop_psrc(gf, gt, xg, job_tag, inv_type, inv_acc, *,
        idx, sfw, sfw_hvp = None, path_sp, psel, fsel, fselc, eig, finished_tags):
    tag = f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_psrc: idx={idx} xg={xg}", job_tag, inv_type, inv_acc)
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig)
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    src = q.mk_point_src(geo, xg)
    fn_sp = os.path.join(path_sp, f"{tag}.lat")
    prop = compute_prop(inv, src, tag = tag, sfw = sfw, fn_sp = fn_sp, psel = psel, fsel = fsel, fselc = fselc)
    if sfw_hvp is not None:
        chvp_16 = q.contract_chvp_16(prop, prop)
        chvp_16.save_float_from_double(sfw_hvp, tag)
    prop_gt = gt * prop
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    prop_gt.glb_sum_tslice().save(get_save_path(fn_spw))

@q.timer_verbose
def compute_prop_psrc_all(gf, gt, pi, job_tag, inv_type, *,
        path_s, path_hvp = None, path_sp, psel, fsel, fselc, eig):
    if q.does_file_exist_sync_node(get_save_path(path_s + ".acc.partial")):
        q.qrename_info(get_save_path(path_s + ".acc.partial"), get_save_path(path_s + ".acc"))
    finished_tags = q.properly_truncate_fields_sync_node(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 8, ])
    sfw_hvp = None
    if path_hvp is not None:
        finished_hvp_tags = q.properly_truncate_fields_sync_node(get_save_path(path_hvp + ".acc"))
        assert finished_tags == finished_hvp_tags
        sfw_hvp = q.open_fields(get_save_path(path_hvp + ".acc"), "a", [ 1, 1, 1, 4, ])
    for inv_acc in [ 2, 1, 0 ]:
        for p in pi:
            idx, xg, inv_type_p, inv_acc_p = p
            if inv_type_p == inv_type and inv_acc_p == inv_acc:
                compute_prop_psrc(gf, gt, xg, job_tag, inv_type, inv_acc,
                        idx = idx, sfw = sfw, sfw_hvp = sfw_hvp, path_sp = path_sp,
                        psel = psel, fsel = fsel, fselc = fselc, eig = eig,
                        finished_tags = finished_tags)
        q.displayln_info(pprint.pformat(q.list_cache()))
        q.clean_cache(q.cache_inv)
        q.displayln_info(pprint.pformat(q.list_cache()))
    sfw_hvp.close()
    sfw.close()
    q.qrename_info(get_save_path(path_hvp + ".acc"), get_save_path(path_hvp))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

@q.timer_verbose
def run_gf(job_tag, traj):
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    if path_gf is None:
        if job_tag[:5] == "test-":
            gf = ru.mk_sample_gauge_field(job_tag, f"{traj}")
            q.qmkdir_info(get_save_path(f"configs"))
            q.qmkdir_info(get_save_path(f"configs/{job_tag}"))
            path_gf = get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
            # gf.save(path_gf)
            qg.save_gauge_field(gf, path_gf)
        else:
            assert False
    get_gf = ru.load_config_lazy(job_tag, path_gf)
    return get_gf

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

@q.timer_verbose
def run_gt(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    path_gt = get_load_path(f"gauge-transform/{job_tag}/traj={traj}.field")
    if path_gt is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-gauge_fix_coulomb"):
            gf = get_gf()
            q.qmkdir_info(get_save_path(f"gauge-transform"))
            q.qmkdir_info(get_save_path(f"gauge-transform/{job_tag}"))
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

@q.timer_verbose
def run_psel(job_tag, traj):
    path_psel = get_load_path(f"point-selection/{job_tag}/traj={traj}.txt")
    if path_psel is None:
        q.qmkdir_info(get_save_path(f"point-selection"))
        q.qmkdir_info(get_save_path(f"point-selection/{job_tag}"))
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

@q.timer_verbose
def run_fsel(job_tag, traj, get_psel):
    if get_psel is None:
        return None
    path_fsel = get_load_path(f"field-selection/{job_tag}/traj={traj}.field")
    total_site = ru.get_total_site(job_tag)
    n_per_tslice = total_site[0] * total_site[1] * total_site[2] // 16
    if path_fsel is None:
        q.qmkdir_info(get_save_path(f"field-selection"))
        q.qmkdir_info(get_save_path(f"field-selection/{job_tag}"))
        fsel = mk_rand_fsel(job_tag, traj, n_per_tslice)
        fsel.save(get_save_path(f"field-selection/{job_tag}/traj={traj}.field"))
        fselc = mk_fselc(fsel, get_psel())
        return lambda : ( fsel, fselc, )
    else:
        @q.timer_verbose
        def load_fsel():
            fsel = q.FieldSelection()
            fsel.load(path_fsel, n_per_tslice)
            fselc = mk_fselc(fsel, get_psel())
            return fsel, fselc
        return q.lazy_call(load_fsel)

@q.timer
def run_prop_wsrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel):
    if None in [ get_gf, get_eig, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-wsrc-light/{job_tag}/traj={traj}") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-light"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        fsel, fselc = get_fsel()
        q.qmkdir_info(get_save_path(f"wall-src-info-light"))
        q.qmkdir_info(get_save_path(f"wall-src-info-light/{job_tag}"))
        q.qmkdir_info(get_save_path(f"prop-wsrc-light"))
        q.qmkdir_info(get_save_path(f"prop-wsrc-light/{job_tag}"))
        q.qmkdir_info(get_save_path(f"psel-prop-wsrc-light"))
        q.qmkdir_info(get_save_path(f"psel-prop-wsrc-light/{job_tag}"))
        q.qmkdir_info(get_save_path(f"psel-prop-wsrc-light/{job_tag}/traj={traj}"))
        wi_light = mk_rand_wall_src_info(job_tag, traj, inv_type = 0)
        save_wall_src_info(wi_light, f"wall-src-info-light/{job_tag}/traj={traj}.txt");
        compute_prop_wsrc_all(gf, gt, wi_light, job_tag, inv_type = 0,
                path_s = f"prop-wsrc-light/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-wsrc-light/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

@q.timer
def run_prop_psrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel):
    if None in [ get_gf, get_eig, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-psrc-light/{job_tag}/traj={traj}") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-psrc-light"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        fsel, fselc = get_fsel()
        q.qmkdir_info(get_save_path(f"point-src-info"))
        q.qmkdir_info(get_save_path(f"point-src-info/{job_tag}"))
        q.qmkdir_info(get_save_path(f"prop-psrc-light"))
        q.qmkdir_info(get_save_path(f"prop-psrc-light/{job_tag}"))
        q.qmkdir_info(get_save_path(f"psel-prop-psrc-light"))
        q.qmkdir_info(get_save_path(f"psel-prop-psrc-light/{job_tag}"))
        q.qmkdir_info(get_save_path(f"psel-prop-psrc-light/{job_tag}/traj={traj}"))
        q.qmkdir_info(get_save_path(f"hvp-psrc-light/"))
        q.qmkdir_info(get_save_path(f"hvp-psrc-light/{job_tag}"))
        pi = mk_rand_point_src_info(job_tag, traj, get_psel())
        save_point_src_info(pi, f"point-src-info/{job_tag}/traj={traj}.txt");
        compute_prop_psrc_all(gf, gt, pi, job_tag, inv_type = 0,
                path_s = f"prop-psrc-light/{job_tag}/traj={traj}",
                path_hvp = f"hvp-psrc-light/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-psrc-light/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

@q.timer
def run_prop_wsrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-wsrc-strange/{job_tag}/traj={traj}") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-strange"):
        gf = get_gf()
        gt = get_gt()
        fsel, fselc = get_fsel()
        q.qmkdir_info(get_save_path(f"wall-src-info-strange"))
        q.qmkdir_info(get_save_path(f"wall-src-info-strange/{job_tag}"))
        q.qmkdir_info(get_save_path(f"prop-wsrc-strange"))
        q.qmkdir_info(get_save_path(f"prop-wsrc-strange/{job_tag}"))
        q.qmkdir_info(get_save_path(f"psel-prop-wsrc-strange"))
        q.qmkdir_info(get_save_path(f"psel-prop-wsrc-strange/{job_tag}"))
        q.qmkdir_info(get_save_path(f"psel-prop-wsrc-strange/{job_tag}/traj={traj}"))
        wi_strange = mk_rand_wall_src_info(job_tag, traj, inv_type = 1)
        save_wall_src_info(wi_strange, f"wall-src-info-strange/{job_tag}/traj={traj}.txt");
        compute_prop_wsrc_all(gf, gt, wi_strange, job_tag, inv_type = 1,
                path_s = f"prop-wsrc-strange/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-wsrc-strange/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = None)
        q.release_lock()

@q.timer
def run_prop_psrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-psrc-strange/{job_tag}/traj={traj}") is not None:
        return
    fsel, fselc = get_fsel()
    if q.obtain_lock(f"locks/{job_tag}-{traj}-psrc-strange"):
        gf = get_gf()
        gt = get_gt()
        q.qmkdir_info(get_save_path(f"point-src-info"))
        q.qmkdir_info(get_save_path(f"point-src-info/{job_tag}"))
        q.qmkdir_info(get_save_path(f"prop-psrc-strange"))
        q.qmkdir_info(get_save_path(f"prop-psrc-strange/{job_tag}"))
        q.qmkdir_info(get_save_path(f"psel-prop-psrc-strange"))
        q.qmkdir_info(get_save_path(f"psel-prop-psrc-strange/{job_tag}"))
        q.qmkdir_info(get_save_path(f"psel-prop-psrc-strange/{job_tag}/traj={traj}"))
        q.qmkdir_info(get_save_path(f"hvp-psrc-strange/"))
        q.qmkdir_info(get_save_path(f"hvp-psrc-strange/{job_tag}"))
        pi = mk_rand_point_src_info(job_tag, traj, get_psel())
        save_point_src_info(pi, f"point-src-info/{job_tag}/traj={traj}.txt");
        compute_prop_psrc_all(gf, gt, pi, job_tag, inv_type = 1,
                path_s = f"prop-psrc-strange/{job_tag}/traj={traj}",
                path_hvp = f"hvp-psrc-strange/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-psrc-strange/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = None)
        q.release_lock()

@q.timer
def run_job(job_tag, traj):
    if check_job(job_tag, traj):
        return
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    get_eig = run_eig(job_tag, traj_gf, get_gf)
    #
    get_psel = run_psel(job_tag, traj)
    get_fsel = run_fsel(job_tag, traj, get_psel)
    assert get_psel is not None
    assert get_fsel is not None
    #
    run_prop_wsrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel)
    run_prop_psrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel)
    run_prop_wsrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel)
    run_prop_psrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel)
    #
    q.timer_display()

rup.dict_params["test-4nt8"]["n_points"] = [
        [ 6, 2, 1, ],
        [ 3, 2, 1, ],
        ]

rup.dict_params["test-4nt16"]["n_points"] = [
        [ 32, 4, 2, ],
        [ 16, 4, 2, ],
        ]

rup.dict_params["48I"]["n_points"] = [
        [ 2048, 64, 16, ],
        [ 1024, 64, 16, ],
        ]

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["48I"]["trajs"] = list(range(500, 3000, 5))

rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10

rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8",
        "test-4nt16",
        # "test-8nt16",
        # "test-16nt32",
        # "test-32nt64",
        # "test-48nt96",
        # "test-64nt128",
        # "test-96nt192",
        # "test-128nt256",
        # "24D",
        # "48I",
        ]

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        q.displayln_info(pprint.pformat(q.list_cache()))
        run_job(job_tag, traj)

qg.end_with_gpt()
