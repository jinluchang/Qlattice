import rbc_ukqcd_params as rup

job_tag = "32Dfine"

dict_params = {}

rup.dict_params[job_tag] = dict_params

dict_params["job_tag"] = job_tag

dict_params["total_site"] = [ 32, 32, 32, 64, ]

dict_params["load_config_params"] = { "twist_boundary_at_boundary":[ 0.0, 0.0, 0.0, -0.5, ] }

def mk_fermion_params(inv_type, inv_acc):
    params = {}
    params["M5"] = 1.8
    params["boundary_phases"] = [1.0, 1.0, 1.0, 1.0] # twist boundary after loading configuration
    params["b"] = 1.0
    params["c"] = 0.0
    if inv_type == 0:
        params["mass"] = 0.0001
    elif inv_type == 1:
        params["mass"] = 0.045
    else:
        assert False
    params["b"] = 0.5 + 32/12/2
    params["c"] = -0.5 + 32/12/2
    params["Ls"] = 12
    return params

def mk_dict_fermion_params():
    params = {}
    for inv_type in [0, 1,]:
        params[inv_type] = {}
        for inv_acc in [0, 1, 2,]:
            params[inv_type][inv_acc] = mk_fermion_params(inv_type, inv_acc)
    return params

dict_params["fermion_params"] = mk_dict_fermion_params()

def mk_lanc_params(inv_type, inv_acc):
    assert inv_type == 0
    assert inv_acc == 0
    fermion_params = mk_dict_fermion_params()[inv_type][inv_acc]
    pit_params = { "eps": 0.01, "maxiter": 500, "real": True, }
    cheby_params = { "low": 7e-5, "high": 5.5, "order": 400, }
    irl_params = {
            "Nstop": 250,
            "Nk": 270,
            "Nm": 350,
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

dict_params["lanc_params"] = { 0: { 0: mk_lanc_params(0, 0), }, }

def mk_clanc_params(inv_type, inv_acc):
    assert inv_type == 0
    assert inv_acc == 0
    block = [ 2, 2, 2, 2, ]
    nbasis = 250
    cheby_params = { "low": 0.0018, "high": 5.5, "order": 100, }
    irl_params = {
            "Nstop": 2000,
            "Nk": 2100,
            "Nm": 2600,
            "resid": 1e-8,
            "betastp": 0.0,
            "maxiter": 20,
            "Nminres": 0,
            # "maxapply": 100
            }
    smoother_params = { "eps": 1e-6, "maxiter": 25, }
    save_params = { "nsingle": 100, "mpi": [ 2, 2, 2, 4, ], }
    return {
            "block": block,
            "nbasis": nbasis,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            "smoother_params": smoother_params,
            "save_params": save_params,
            }

dict_params["clanc_params"] = { 0: { 0: mk_clanc_params(0, 0), }, }
