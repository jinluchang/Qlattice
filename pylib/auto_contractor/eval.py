import qlat_gpt as qg
import qlat as q
import gpt as g
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import numpy as np
import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [ "results" ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer
def mk_sample_gauge_field(job_tag, fn):
    rs = q.RngState(f"seed {job_tag} {fn}").split("mk_sample_gauge_field")
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    gf = q.GaugeField(geo)
    gf.set_rand(rs, sigma = 0.25, n_step = 12)
    for i in range(7):
        q.gf_wilson_flow_step(gf, 0.05)
    gf.unitarize()
    return gf

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
    ru.save_ceig(get_save_path(path), eig, job_tag, inv_type, inv_acc);
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
    for inv_acc in [0, 1,]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig, mpi_split = False, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} with eig")
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, mpi_split = False, timer = False) * src
        sol -= sol_ref
        q.displayln_info(f"sol diff norm {sol.qnorm()} inv_acc={inv_acc} without eig")

@q.timer
def compute_prop(inv, src, *, tag, sfw):
    sol = inv * src
    sol.save_double(sfw, tag)

@q.timer
def compute_prop_psrc(gf, xg, job_tag, inv_type, inv_acc, *, idx, sfw, eig, finished_tags):
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
    compute_prop(inv, src, tag = tag, sfw = sfw)

def get_all_points(total_site, *, tslice = None):
    all_points = []
    if tslice is None:
        for t in range(total_site[3]):
            all_points += get_all_points(total_site, tslice = t)
        return all_points
    t = tslice
    for x in range(total_site[0]):
        for y in range(total_site[1]):
            for z in range(total_site[2]):
                all_points.append([x, y, z, t,])
    return all_points

@q.timer
def compute_prop_psrc_all(gf, job_tag, inv_type, *, path_s, eig):
    finished_tags = q.properly_truncate_fields_sync_node(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 2 ])
    total_site = ru.get_total_site(job_tag)
    inv_acc = 2
    for idx, xg in enumerate(get_all_points(total_site)):
        compute_prop_psrc(gf, xg, job_tag, inv_type, inv_acc,
                idx = idx, sfw = sfw, eig = eig,
                finished_tags = finished_tags)
    q.clean_cache(q.cache_inv)
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

prop_cache = q.mk_cache("prop_cache") # prop_cache[flavor][src_p] = prop

@q.timer
def get_prop_psrc(flavor : str, xg):
    # call load_prop_psrc_all(flavor, path_s) first
    return prop_cache[flavor][f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]})"]

@q.timer
def get_prop_psrc_psnk(flavor : str, xg_src, xg_snk):
    return np.array(get_prop_psrc(flavor, xg_src).get_elem(xg_snk)).reshape((4, 4, 3, 3,))

@q.timer
def load_prop_psrc_all(total_site, flavor : str, path_s : str):
    cache = q.mk_cache(flavor, cache = prop_cache)
    get_all_points(total_site)
    if flavor in ["l", "u", "d",]:
        inv_type = 0
    elif flavor in ["s",]:
        inv_type = 1
    else:
        assert False
    inv_acc = 2
    sfr = q.open_fields(get_load_path(path_s), "r")
    for idx, xg in enumerate(get_all_points(total_site)):
        q.displayln_info(f"load_prop_psrc_all: idx={idx} xg={xg} flavor={flavor} path_s={path_s}")
        prop = q.Prop()
        tag = f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc}"
        prop.load_double(sfr, tag)
        # convert to GPT/Grid mspincolor wilson matrix order
        prop_msc = q.Prop()
        q.convert_mspincolor_from_wm_prop(prop_msc, prop)
        cache[f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]})"] = prop_msc
    sfr.close()

@q.timer
def run_job(job_tag, traj):
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(rup.dict_params[job_tag])
    #
    q.qmkdir_info(f"locks")
    q.qmkdir_info(get_save_path(f""))
    #
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    q.displayln_info("geo.show() =", geo.show())
    #
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    if path_gf is None:
        if job_tag[:5] == "test-":
            gf = ru.mk_sample_gauge_field(job_tag, f"{traj}")
            q.qmkdir_info(get_save_path(f"configs"))
            q.qmkdir_info(get_save_path(f"configs/{job_tag}"))
            path_gf = get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
            gf.save(path_gf)
        else:
            assert False
    gf = ru.load_config(job_tag, path_gf)
    gf.show_info()
    #
    get_eig = None
    if q.obtain_lock(f"locks/{job_tag}-{traj}-compute-eig"):
        q.qmkdir_info(get_save_path(f"eig"))
        q.qmkdir_info(get_save_path(f"eig/{job_tag}"))
        get_eig = compute_eig(gf, job_tag, inv_type = 0, path = f"eig/{job_tag}/traj={traj}")
        q.release_lock()
    #
    for inv_type in [0, 1,]:
        if get_load_path(f"prop-psrc-{inv_type}/{job_tag}/traj={traj}") is None:
            if inv_type == 0 and get_eig is None:
                continue
            if q.obtain_lock(f"locks/{job_tag}-{traj}-compute-prop-psrc-all-{inv_type}"):
                q.qmkdir_info(get_save_path(f"prop-psrc-{inv_type}"))
                q.qmkdir_info(get_save_path(f"prop-psrc-{inv_type}/{job_tag}"))
                if inv_type == 0:
                    eig = get_eig()
                else:
                    eig = None
                compute_prop_psrc_all(gf, job_tag, inv_type, path_s = f"prop-psrc-{inv_type}/{job_tag}/traj={traj}", eig = eig)
                q.release_lock()
    #
    path_prop_list = [ get_load_path(f"prop-psrc-{inv_type}/{job_tag}/traj={traj}") for inv_type in [0, 1,] ]
    if all(map(lambda x : x is not None, path_prop_list)):
        load_prop_psrc_all(total_site, "l", f"prop-psrc-0/{job_tag}/traj={traj}")
        load_prop_psrc_all(total_site, "s", f"prop-psrc-1/{job_tag}/traj={traj}")
        q.displayln_info(get_prop_psrc_psnk("l", [1, 2, 3, 2], [1, 2, 3, 2]))
    #
    q.clean_cache()

if __name__ == "__main__":
    qg.begin_with_gpt()
    job_tag = "test-4nt16"
    traj = 1000
    rup.dict_params[job_tag]["fermion_params"][0][2] = rup.dict_params[job_tag]["fermion_params"][0][0]
    rup.dict_params[job_tag]["fermion_params"][1][2] = rup.dict_params[job_tag]["fermion_params"][1][0]
    rup.dict_params[job_tag]["load_config_params"]["twist_boundary_at_boundary"] = [0.0, 0.0, 0.0, -0.5,]
    run_job(job_tag, traj)
    q.timer_display()
    qg.end_with_gpt()
