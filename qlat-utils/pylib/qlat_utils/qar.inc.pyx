def get_qar_multi_vol_max_size(size = None):
    """
    Parameter controls the size of a single `qar` file in number of bytes. Note, `qar` never splits a single file into multiple `qar` volume.
    """
    cdef long* p_size
    if size is not None:
        assert isinstance(size, int)
        p_size = &cc.cc_get_qar_multi_vol_max_size()
        p_size[0] = <long>size
    return cc.cc_get_qar_multi_vol_max_size()

def does_regular_file_exist_qar(const cc.std_string& path):
    return cc.cc_does_regular_file_exist_qar(path)

def does_file_exist_qar(const cc.std_string& path):
    return cc.cc_does_file_exist_qar(path)

def qcat(const cc.std_string& path):
    """Return contents of file as `str`"""
    return <str>cc.cc_qcat(path)

def qcat_bytes(const cc.std_string& path):
    """Return contents of file as `bytes`"""
    return <bytes>cc.cc_qcat(path)

@timer
def qar_create(const cc.std_string& path_qar, const cc.std_string& path_folder,
               *, const cc.bool is_remove_folder_after = False):
    return cc.cc_qar_create(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract(const cc.std_string& path_qar, const cc.std_string& path_folder,
                *, const cc.bool is_remove_qar_after = False):
    return cc.cc_qar_extract(path_qar, path_folder, is_remove_qar_after)

@timer
def qcopy_file(const cc.std_string& path_src, const cc.std_string& path_dst):
    return cc.cc_qcopy_file(path_src, path_dst)

@timer
def list_qar(const cc.std_string& path_qar):
    cdef list l = cc.cc_list_qar(path_qar)
    return [ str(fn) for fn in l ]
