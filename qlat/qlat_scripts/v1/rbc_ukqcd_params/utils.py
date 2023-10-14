from .. import rbc_ukqcd_params as rup

from . import dict_params

import qlat as q

def get_param(*keys, default = None, dict_params = None):
    if dict_params is None:
        dict_params = rup.dict_params
    d = dict_params
    for key in keys:
        if key in d:
            d = d[key]
        else:
            return default
    return d

def set_param(*keys, value, dict_params = None):
    if dict_params is None:
        dict_params = rup.dict_params
    assert len(keys) >= 1
    d = dict_params
    for key in keys[:-1]:
        if key not in d:
            d[key] = dict()
        d = d[key]
    d[keys[-1]] = value

def get_total_site(job_tag : str):
    return get_param(job_tag, "total_site")

@q.timer_verbose
def mk_sample_gauge_field(job_tag, fn):
    rs = q.RngState(f"seed {job_tag} {fn}").split("mk_sample_gauge_field")
    total_site = get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    gf = q.GaugeField(geo)
    gf.set_rand(rs, sigma = 0.25, n_step = 4)
    for i in range(4):
        q.gf_wilson_flow_step(gf, 0.05)
    gf.unitarize()
    return gf

@q.timer_verbose
def mk_sample_gauge_field_v2(total_site, tag):
    rs = q.RngState(f"seed {total_site} {tag}").split("mk_sample_gauge_field")
    geo = q.Geometry(total_site, 1)
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split("gf.set_rand"), sigma = 0.25, n_step = 16)
    gf.unitarize()
    beta = 6.0
    ga = q.GaugeAction(beta)
    for traj in range(500):
        q.run_hmc_pure_gauge(gf, ga, traj, rs.split("run_hmc_pure_gauge"),
                             n_step = 6, is_always_accept = True)
    gf.unitarize()
    return gf

@q.timer_verbose
def mk_sample_gauge_field_v3(job_tag, fn):
    """
    depends on
    total_site = get_param(job_tag, "total_site")
    rand_n_step = get_param(job_tag, "mk_sample_gauge_field", "rand_n_step", default = 16)
    rand_sigma = get_param(job_tag, "mk_sample_gauge_field", "rand_sigma", default = 0.25)
    flow_n_step = get_param(job_tag, "mk_sample_gauge_field", "flow_n_step", default = 4)
    flow_size = get_param(job_tag, "mk_sample_gauge_field", "flow_size", default = 0.05)
    hmc_n_traj = get_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj", default = 500)
    hmc_beta = get_param(job_tag, "mk_sample_gauge_field", "hmc_beta", default = 6.0)
    hmc_n_step = get_param(job_tag, "mk_sample_gauge_field", "hmc_n_step", default = 6)
    hmc_is_always_accept = get_param(job_tag, "mk_sample_gauge_field", "hmc_is_always_accept", default = True)
    """
    rs = q.RngState(f"seed {job_tag} {fn}").split("mk_sample_gauge_field")
    total_site = get_param(job_tag, "total_site")
    rand_n_step = get_param(job_tag, "mk_sample_gauge_field", "rand_n_step", default = 16)
    rand_sigma = get_param(job_tag, "mk_sample_gauge_field", "rand_sigma", default = 0.25)
    flow_n_step = get_param(job_tag, "mk_sample_gauge_field", "flow_n_step", default = 4)
    flow_size = get_param(job_tag, "mk_sample_gauge_field", "flow_size", default = 0.05)
    hmc_n_traj = get_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj", default = 500)
    hmc_beta = get_param(job_tag, "mk_sample_gauge_field", "hmc_beta", default = 6.0)
    hmc_n_step = get_param(job_tag, "mk_sample_gauge_field", "hmc_n_step", default = 6)
    hmc_is_always_accept = get_param(job_tag, "mk_sample_gauge_field", "hmc_is_always_accept", default = True)
    geo = q.Geometry(total_site, 1)
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split("gf.set_rand"), sigma = rand_sigma, n_step = rand_n_step)
    gf.unitarize()
    for i in range(flow_n_step):
        q.gf_wilson_flow_step(gf, flow_size)
    gf.unitarize()
    ga = q.GaugeAction(hmc_beta)
    for traj in range(hmc_n_traj):
        q.run_hmc_pure_gauge(gf, ga, traj, rs.split("run_hmc_pure_gauge"),
                             n_step = hmc_n_step,
                             is_always_accept = hmc_is_always_accept)
    gf.unitarize()
    return gf

@q.timer_verbose
def load_config(job_tag : str, fn : str):
    if not q.does_file_exist_qar_sync_node(fn):
        raise Exception(f"load_config '{fn}' does not exist.")
    gf = q.GaugeField()
    gf.load(fn)
    # gf = qg.load_gauge_field(fn)
    if (job_tag in dict_params) and ("load_config_params" in dict_params[job_tag]):
        params = dict_params[job_tag]["load_config_params"]
        twist_boundary_at_boundary = params["twist_boundary_at_boundary"]
        for mu in range(4):
            lmom = twist_boundary_at_boundary[mu]
            if lmom != 0.0:
                q.displayln_info(f"load_config fn='{fn}' twist_boundary_at_boundary lmom={lmom} mu={mu}")
                gf.twist_boundary_at_boundary(lmom, mu)
    gf.show_info()
    return gf

@q.timer_verbose
def load_config_lazy(job_tag : str, fn : str):
    if not q.does_file_exist_qar_sync_node(fn):
        return None
    return q.lazy_call(load_config, job_tag, fn)
