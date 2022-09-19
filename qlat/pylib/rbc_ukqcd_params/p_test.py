import rbc_ukqcd_params as rup

def mk_test_l_t_list():
    lt_list = []
    for l in [ 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, ]:
        for t in [ 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, ]:
            lt_list.append([l, t,])
    return lt_list

def mk_fermion_params(inv_type, inv_acc):
    params = {}
    params["M5"] = 1.8
    params["boundary_phases"] = [ 1.0, 1.0, 1.0, 1.0, ] # twist boundary after loading configuration
    if inv_acc in [ 0, 1, 2, ]:
        params["b"] = 1.5
        params["c"] = 0.5
        params["Ls"] = 8
    else:
        assert False
    if inv_type == 0:
        params["mass"] = 0.01
    elif inv_type == 1:
        params["mass"] = 0.04
    elif inv_type == 2:
        params["mass"] = 0.20
    else:
        assert False
    return params

def mk_dict_fermion_params():
    params = {}
    for inv_type in [ 0, 1, 2, ]:
        params[inv_type] = {}
        for inv_acc in [ 0, 1, 2, ]:
            params[inv_type][inv_acc] = mk_fermion_params(inv_type, inv_acc)
    return params

def mk_lanc_params(job_tag, inv_type, inv_acc):
    assert inv_type in [ 0, 1, ]
    assert inv_acc == 0
    fermion_params = mk_dict_fermion_params()[inv_type][inv_acc]
    pit_params = { "eps": 0.01, "maxiter": 500, "real": True }
    if job_tag == "test-4nt8":
        cheby_params = {"low": 0.10, "high": 5.5, "order": 50}
    elif job_tag == "test-4nt16":
        cheby_params = {"low": 0.05, "high": 5.5, "order": 50}
    else:
        cheby_params = {"low": 0.05, "high": 5.5, "order": 50}
    irl_params = {
            "Nstop": 20,
            "Nk": 25,
            "Nm": 30,
            "resid": 1e-8,
            "betastp": 0.0,
            "maxiter": 20,
            "Nminres": 1,
            # "maxapply": 100
            }
    return {
            "fermion_params": fermion_params,
            "pit_params": pit_params,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            }

def mk_clanc_params(job_tag, inv_type, inv_acc):
    assert inv_type in [ 0, 1, ]
    assert inv_acc == 0
    block = [ 2, 2, 2, 2 ]
    nbasis = 20
    if inv_type == 0:
        cheby_params = { "low": 0.20, "high": 5.5, "order": 50, }
    elif inv_type == 1:
        cheby_params = { "low": 0.30, "high": 5.5, "order": 50, }
    else:
        assert False
    irl_params = {
            "Nstop": 30,
            "Nk": 35,
            "Nm": 40,
            "resid": 1e-8,
            "betastp": 0.0,
            "maxiter": 20,
            "Nminres": 1,
            # "maxapply": 100
            }
    smoother_params = { "eps": 1e-8, "maxiter": 10, }
    save_params = { "nsingle": 10, "mpi": [ 1, 1, 1, 4, ], }
    return {
            "block": block,
            "nbasis": nbasis,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            "smoother_params": smoother_params,
            "save_params": save_params,
            }

def setup_params():
    for l, t in mk_test_l_t_list():
        job_tag = f"test-{l}nt{t}"
        dict_params = {}
        rup.dict_params[job_tag] = dict_params
        dict_params["job_tag"] = job_tag
        dict_params["total_site"] = [ l, l, l, t, ]
        dict_params["load_config_params"] = { "twist_boundary_at_boundary" : [ 0.0, 0.0, 0.0, -0.5, ], }
        dict_params["fermion_params"] = mk_dict_fermion_params()
        if job_tag == "test-4nt8":
            dict_lanc_params = { 0: { 0: mk_lanc_params(job_tag, 0, 0), }, }
            dict_clanc_params = { 0: { 0: mk_clanc_params(job_tag, 0, 0), }, }
        else:
            dict_lanc_params = { 0: { 0: mk_lanc_params(job_tag, 0, 0), }, 1: { 0: mk_lanc_params(job_tag, 1, 0), }, }
            dict_clanc_params = { 0: { 0: mk_clanc_params(job_tag, 0, 0), }, 1: { 0: mk_clanc_params(job_tag, 1, 0), }, }
        dict_params["lanc_params"] = dict_lanc_params
        dict_params["clanc_params"] = dict_clanc_params

setup_params()
