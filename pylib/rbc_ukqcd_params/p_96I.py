import rbc_ukqcd_params as rup

job_tag = "96I"

dict_params = {}

rup.dict_params[job_tag] = dict_params

dict_params["job_tag"] = job_tag

dict_params["total_site"] = [ 96, 96, 96, 192, ]

dict_params["load_config_params"] = { "twist_boundary_at_boundary":[ 0.0, 0.0, 0.0, -0.5, ] }
