import rbc_ukqcd_params as rup

import qlat as q

def get_total_site(job_tag : str):
    return rup.dict_params[job_tag]["total_site"]

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
def load_config(job_tag : str, fn : str):
    if not q.does_file_exist_qar_sync_node(fn):
        raise Exception(f"load_config '{fn}' does not exist.")
    gf = q.GaugeField()
    gf.load(fn)
    # gf = qg.load_gauge_field(fn)
    if (job_tag in rup.dict_params) and ("load_config_params" in rup.dict_params[job_tag]):
        params = rup.dict_params[job_tag]["load_config_params"]
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
