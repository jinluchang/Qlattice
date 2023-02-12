def flush():
    cc.flush()

def qtouch(const cc.std_string& path, content = None):
    if content is None:
        return cc.cc_qtouch(path)
    else:
        return cc.cc_qtouch(path, content)

def qtouch_info(const cc.std_string& path, content = None):
    if content is None:
        return cc.cc_qtouch_info(path)
    else:
        return cc.cc_qtouch_info(path, content)

def qappend(const cc.std_string& path, const cc.std_string& content):
    return cc.cc_qappend(path, content)

def qappend_info(const cc.std_string& path, const cc.std_string& content):
    return cc.cc_qappend_info(path, content)

def qrename(const cc.std_string& old_path, const cc.std_string& new_path):
    return cc.cc_qrename(old_path, new_path)

def qrename_info(const cc.std_string& old_path, const cc.std_string& new_path):
    return cc.cc_qrename_info(old_path, new_path)

def qremove(const cc.std_string& path):
    return cc.cc_qremove(path)

def qremove_all(const cc.std_string& path):
    return cc.cc_qremove_all(path)

def qmkdir(const cc.std_string& path):
    return cc.cc_qmkdir(path)

def qmkdir_info(const cc.std_string& path):
    return cc.cc_qmkdir_info(path)

def qmkdir_p(const cc.std_string& path):
    return cc.cc_qmkdir_p(path)

def qmkdir_p_info(const cc.std_string& path):
    return cc.cc_qmkdir_p_info(path)

def does_file_exist(const cc.std_string& path):
    return cc.cc_does_file_exist(path)

def is_directory(const cc.std_string& path):
    return cc.cc_is_directory(path)

def is_regular_file(const cc.std_string& path):
    return cc.cc_is_regular_file(path)

def qls(const cc.std_string& path,
        const cc.bool is_sort = True):
    cdef list l = cc.cc_qls(path, is_sort)
    return [ str(v) for v in l ]

def qls_all(const cc.std_string& path,
            const cc.bool is_folder_before_files = False,
            const cc.bool is_sort = True):
    cdef list l = cc.cc_qls_all(path, is_folder_before_files, is_sort)
    return [ str(v) for v in l ]

def compute_crc32(const cc.std_string& path):
    return cc.cc_compute_crc32(path)

def qload_datatable(const cc.std_string& path, const cc.bool is_par = False):
    return cc.cc_qload_datatable(path, is_par)

def check_all_files_crc32_info(const cc.std_string& path):
    return cc.cc_check_all_files_crc32_info(path)

