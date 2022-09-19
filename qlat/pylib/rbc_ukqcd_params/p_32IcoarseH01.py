import rbc_ukqcd_params as rup

# Christoph ensemble 1, ml = 0.0025, ms = 0.0362, Ls = 24, beta = 2.13

job_tag = "32IcoarseH01"

dict_params = {}

rup.dict_params[job_tag] = dict_params

dict_params["job_tag"] = job_tag

dict_params["total_site"] = [ 32, 32, 32, 64, ]

dict_params["load_config_params"] = { "twist_boundary_at_boundary": [ 0.0, 0.0, 0.0, -0.5, ] }

def mk_fermion_params(inv_type, inv_acc):
    params = {}
    params["M5"] = 1.8
    params["boundary_phases"] = [ 1.0, 1.0, 1.0, 1.0, ] # twist boundary after loading configuration
    params["b"] = 1.5
    params["c"] = 0.5
    if inv_type == 0:
        params["mass"] = 0.0025
    elif inv_type == 1:
        params["mass"] = 0.0362
    elif inv_type == 2:
        params["mass"] = 0.2
    else:
        assert False
    if inv_acc == 0 or inv_acc == 1:
        params["b"] = 1.0
        params["c"] = 0.0
        params["omega"] = [
                1.4789834351796358,
                1.347049274947458,
                1.1273467425761714,
                0.891777252638092,
                0.6798157283073448,
                0.506488728896523,
                0.371660686773301,
                0.26925807172929905,
                0.19187135204541664,
                0.1516448422854592,
                0.13000288981245012 + 0.04732386585347784j,
                0.13000288981245012 - 0.04732386585347784j,
                0.08653188911012571 + 0.10514692271165546j,
                0.08653188911012571 - 0.10514692271165546j,
                ]
    elif inv_acc == 2:
        params["Ls"] = 24
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
