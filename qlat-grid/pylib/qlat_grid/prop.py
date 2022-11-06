import cqlat_grid as cgi
import qlat as q

def save_prop_float(prop, path):
    cgi.save_prop_float(prop, path + ".partial")
    q.qrename_info(path + ".partial", path)

def load_prop_float(prop, path):
    # prop need to have the right geo as input
    return cgi.load_prop_float(prop, path)
