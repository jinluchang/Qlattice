import rbc_ukqcd_params as rup

# Christoph ensemble 9, ml = 0.00372, ms = 0.0257, Ls = 12, beta = 2.25

job_tag = "32IH01"

dict_params = {}

rup.dict_params[job_tag] = dict_params

dict_params["job_tag"] = job_tag

dict_params["total_site"] = [ 32, 32, 32, 64, ]

dict_params["load_config_params"] = { "twist_boundary_at_boundary":[ 0.0, 0.0, 0.0, -0.5, ] }

def mk_fermion_params(inv_type, inv_acc):
    params = {}
    params["M5"] = 1.8
    params["boundary_phases"] = [ 1.0, 1.0, 1.0, 1.0, ] # twist boundary after loading configuration
    params["b"] = 1.5
    params["c"] = 0.5
    if inv_type == 0:
        params["mass"] = 0.00372
    elif inv_type == 1:
        params["mass"] = 0.0257
    elif inv_type == 2:
        params["mass"] = 0.31
    else:
        assert False
    params["Ls"] = 12
    return params

def mk_dict_fermion_params():
    params = {}
    for inv_type in [ 0, 1, 2, ]:
        params[inv_type] = {}
        for inv_acc in [ 0, 1, 2, ]:
            params[inv_type][inv_acc] = mk_fermion_params(inv_type, inv_acc)
    return params

dict_params["fermion_params"] = mk_dict_fermion_params()
