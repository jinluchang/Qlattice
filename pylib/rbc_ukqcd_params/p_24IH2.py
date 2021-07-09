import rbc_ukqcd_params as rup

job_tag = "24IH2"

rup.dict_total_site[job_tag] = [ 24, 24, 24, 64, ]

load_config_params = { "twist_boundary_at_boundary":[ 0.0, 0.0, 0.0, -0.5, ] }

rup.dict_load_config_params[job_tag] = load_config_params

def mk_fermion_params(inv_type, inv_acc):
    params = {}
    params["M5"] = 1.8
    params["boundary_phases"] = [1.0, 1.0, 1.0, 1.0] # twist boundary after loading configuration
    params["b"] = 1.0
    params["c"] = 0.0
    params["Ls"] = 16
    if inv_type == 0:
        params["mass"] = 0.01
    elif inv_type == 1:
        params["mass"] = 0.04
    else:
        assert False
    return params

def mk_dict_fermion_params():
    params = {}
    params[inv_type][inv_acc] = {}
    for inv_type in [0, 1,]:
        params[inv_type] = {}
        for inv_acc in [0, 1, 2,]:
            params[inv_type][inv_acc] = mk_fermion_params(inv_type, inv_acc)
    return params

dict_fermion_params = mk_dict_fermion_params()

rup.dict_fermion_params[job_tag] = dict_fermion_params

def mk_lanc_params(inv_type, inv_acc):
    assert inv_type == 0
    assert inv_acc == 0
    fermion_params = dict_fermion_params[inv_type][inv_acc]
    pit_params = { "eps": 0.01, "maxiter": 500, "real": True }
    cheby_params = {"low": 0.001, "high": 5.5, "order": 100}
    irl_params = {
            "Nstop": 250,
            "Nk": 260,
            "Nm": 300,
            "resid": 1e-8,
            "betastp": 0.0,
            "maxiter": 20,
            "Nminres": 3,
            # "maxapply": 100
            }
    return {
            "fermion_params": fermion_params,
            "pit_params": pit_params,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            }

def mk_clanc_params(inv_type, inv_acc):
    assert inv_type == 0
    assert inv_acc == 0
    block = [ 2, 2, 2, 2, ]
    nbasis = 250
    cheby_params = {"low": 0.0025, "high": 5.5, "order": 100}
    irl_params = {
            "Nstop": 500,
            "Nk": 510,
            "Nm": 550,
            "resid": 1e-8,
            "betastp": 0.0,
            "maxiter": 20,
            "Nminres": 4,
            #    "maxapply" : 100
            }
    smoother_params = {"eps": 1e-6, "maxiter": 10}
    save_params = {"nsingle": 100, "mpi": [ 1, 1, 1, 8, ]}
    return {
            "block": block,
            "nbasis": nbasis,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            "smoother_params": smoother_params,
            "save_params": save_params,
            }

dict_lanc_params = { 0:{ 0:mk_lanc_params(0, 0) } }
dict_clanc_params = { 0:{ 0:mk_clanc_params(0, 0) } }

rup.dict_lanc_params[job_tag] = dict_lanc_params
rup.dict_clanc_params[job_tag] = dict_clanc_params
