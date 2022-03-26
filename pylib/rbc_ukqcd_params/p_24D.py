import rbc_ukqcd_params as rup

job_tag = "24D"

dict_params = {}

rup.dict_params[job_tag] = dict_params

dict_params["job_tag"] = job_tag

dict_params["total_site"] = [ 24, 24, 24, 64, ]

dict_params["load_config_params"] = { "twist_boundary_at_boundary":[ 0.0, 0.0, 0.0, -0.5, ] }

def mk_fermion_params(inv_type, inv_acc):
    params = {}
    params["M5"] = 1.8
    params["boundary_phases"] = [ 1.0, 1.0, 1.0, 1.0, ] # twist boundary after loading configuration
    params["b"] = 2.5
    params["c"] = 1.5
    if inv_type == 0:
        params["mass"] = 0.00107
    elif inv_type == 1:
        params["mass"] = 0.0850
    elif inv_type == 2:
        params["mass"] = 0.25
    else:
        assert False
    if inv_acc == 0 or inv_acc == 1:
        params["b"] = 1.0
        params["c"] = 0.0
        params["omega"] = [
                1.0903256131299373,
                0.9570283702230611,
                0.7048886040934104,
                0.48979921782791747,
                0.328608311201356,
                0.21664245377015995,
                0.14121112711957107,
                0.0907785101745156,
                0.05608303440064219 - 0.007537158177840385j,
                0.05608303440064219 + 0.007537158177840385j,
                0.0365221637144842 - 0.03343945161367745j,
                0.0365221637144842 + 0.03343945161367745j,
                ]
    else:
        params["Ls"] = 24
    return params

def mk_dict_fermion_params():
    params = {}
    for inv_type in [ 0, 1, 2, ]:
        params[inv_type] = {}
        for inv_acc in [ 0, 1, 2, ]:
            params[inv_type][inv_acc] = mk_fermion_params(inv_type, inv_acc)
    return params

dict_params["fermion_params"] = mk_dict_fermion_params()

dict_params[f"cg_params-0-2"] = {}
dict_params[f"cg_params-0-2"]["maxiter"] = 300

dict_params[f"cg_params-1-0"] = {}
dict_params[f"cg_params-1-0"]["maxiter"] = 250

dict_params[f"cg_params-1-1"] = {}
dict_params[f"cg_params-1-1"]["maxiter"] = 250

dict_params[f"cg_params-1-2"] = {}
dict_params[f"cg_params-1-2"]["maxiter"] = 300

def mk_lanc_params(inv_type, inv_acc):
    assert inv_acc == 0
    if inv_type == 0:
        c_low = 2.5e-4
        n_stop = 1000
        n_keep = n_stop + 40
        n_max = n_stop + 200
    elif inv_type == 1:
        c_low = 1.3e-4
        n_stop = 500
        n_keep = n_stop + 40
        n_max = n_stop + 200
    else:
        assert False
    fermion_params = mk_dict_fermion_params()[inv_type][inv_acc]
    pit_params = { "eps": 0.01, "maxiter": 500, "real": True }
    cheby_params = {"low": c_low, "high": 5.5, "order": 200}
    irl_params = {
            "Nstop": n_stop,
            "Nk": n_keep,
            "Nm": n_max,
            "resid": 1e-8,
            "betastp": 0.0,
            "maxiter": 20,
            "Nminres": 0,
            # "maxapply": 100
            }
    return {
            "fermion_params": fermion_params,
            "pit_params": pit_params,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            }

def mk_clanc_params(inv_type, inv_acc):
    assert inv_acc == 0
    if inv_type == 0:
        nbasis = 1000
        c_low = 7e-4
        n_stop = 2000
        n_keep = n_stop + 100
        n_max = n_stop + 600
        n_single = 150
    elif inv_type == 1:
        nbasis = 500
        c_low = 2.5e-4
        n_stop = 1000
        n_keep = n_stop + 50
        n_max = n_stop + 300
        n_single = 100
    else:
        assert False
    block = [ 2, 3, 3, 2, ]
    cheby_params = { "low": c_low, "high": 5.5, "order": 100, }
    irl_params = {
            "Nstop": n_stop,
            "Nk": n_keep,
            "Nm": n_max,
            "resid": 1e-8,
            "betastp": 0.0,
            "maxiter": 20,
            "Nminres": 0,
            #    "maxapply" : 100
            }
    smoother_params = { "eps": 1e-6, "maxiter": 25, }
    save_params = { "nsingle": n_single, "mpi": [ 1, 1, 1, 16, ], }
    return {
            "block": block,
            "nbasis": nbasis,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            "smoother_params": smoother_params,
            "save_params": save_params,
            }

dict_params["lanc_params"] = { 0: { 0: mk_lanc_params(0, 0), }, 1: { 0: mk_lanc_params(1, 0), }, }

dict_params["clanc_params"] = { 0: { 0: mk_clanc_params(0, 0), }, 1: { 0: mk_clanc_params(1, 0), }, }
