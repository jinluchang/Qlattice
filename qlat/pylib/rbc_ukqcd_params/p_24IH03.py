import rbc_ukqcd_params as rup

# Christoph ensemble A, ml = 0.002356, ms = 0.03366, Ls = 8, beta = 2.13

job_tag = "24IH03"

dict_params = {}

rup.dict_params[job_tag] = dict_params

dict_params["job_tag"] = job_tag

dict_params["total_site"] = [ 24, 24, 24, 48, ]

dict_params["load_config_params"] = { "twist_boundary_at_boundary":[ 0.0, 0.0, 0.0, -0.5, ] }

def mk_fermion_params(inv_type, inv_acc):
    params = {}
    params["M5"] = 1.8
    params["boundary_phases"] = [ 1.0, 1.0, 1.0, 1.0, ] # twist boundary after loading configuration
    params["b"] = 1.5
    params["c"] = 0.5
    if inv_type == 0:
        params["mass"] = 0.002356
    elif inv_type == 1:
        params["mass"] = 0.03366
    elif inv_type == 2:
        params["mass"] = 0.2
    else:
        assert False
    if inv_acc == 0 or inv_acc == 1 or inv_acc == 2:
        params["Ls"] = 8
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

dict_params["fermion_params"] = mk_dict_fermion_params()
