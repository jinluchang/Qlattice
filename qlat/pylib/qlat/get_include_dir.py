import os
import qlat.field
import qlat_utils as qu

def get_qlat_include():
    return os.path.join(os.path.dirname(qlat.field.__file__), 'include')

def get_include_list():
    return [ get_qlat_include(), ] + qu.get_include_list()
