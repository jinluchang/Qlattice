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
    if inv_acc == 0 or inv_acc == 1:
        params["b"] = 1.0 + 32/12/2
        params["c"] = 0.0 + 32/12/2
        params["Ls"] = 12
    else:
        params["Ls"] = 32
    return params

def mk_dict_fermion_params():
    params = {}
    for inv_type in [0, 1,]:
        params[inv_type] = {}
        for inv_acc in [0, 1, 2,]:
            params[inv_type][inv_acc] = mk_fermion_params(inv_type, inv_acc)
    return params

dict_params["fermion_params"] = mk_dict_fermion_params()
