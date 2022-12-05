import qlat_grid.c as cg

import qlat as q

def save_prop_float(prop, path):
    cg.save_prop_float(prop, path + ".partial")
    q.qrename_info(path + ".partial", path)

def load_prop_float(prop, path):
    # prop need to have the right geo as input
    return cg.load_prop_float(prop, path)
