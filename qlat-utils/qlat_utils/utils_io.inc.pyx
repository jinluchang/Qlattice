def flush():
    cc.flush()

def qtouch(const cc.std_string& path, content = None):
    if content is None:
        return cc.qtouch(path)
    else:
        return cc.qtouch(path, content)

def qtouch_info(const cc.std_string& path, content = None):
    if content is None:
        return cc.qtouch_info(path)
    else:
        return cc.qtouch_info(path, content)

def qappend(const cc.std_string& path, const cc.std_string& content):
    return cc.qappend(path, content)

def qappend_info(const cc.std_string& path, const cc.std_string& content):
    return cc.qappend_info(path, content)

def qrename(const cc.std_string& old_path, const cc.std_string& new_path):
    return cc.qrename(old_path, new_path)

def qrename_info(const cc.std_string& old_path, const cc.std_string& new_path):
    return cc.qrename_info(old_path, new_path)

def qremove(const cc.std_string& path):
    return cc.qremove(path)

def qremove_all(const cc.std_string& path):
    return cc.qremove_all(path)

def qmkdir(const cc.std_string& path):
    return cc.qmkdir(path)

def qmkdir_info(const cc.std_string& path):
    return cc.qmkdir_info(path)

def qmkdir_p(const cc.std_string& path):
    return cc.qmkdir_p(path)

def qmkdir_p_info(const cc.std_string& path):
    return cc.qmkdir_p_info(path)

def is_directory(const cc.std_string& path):
    return cc.is_directory(path)

def is_regular_file(const cc.std_string& path):
    return cc.is_regular_file(path)

def does_file_exist(const cc.std_string& path):
    return cc.does_file_exist(path)

def clear_is_directory_cache():
    return cc.clear_is_directory_cache()

def is_directory_cache(const cc.std_string& path):
    return cc.is_directory_cache(path)

def is_regular_file_cache(const cc.std_string& path):
    return cc.is_regular_file_cache(path)

def does_file_exist_cache(const cc.std_string& path):
    return cc.does_file_exist_cache(path)

def qls(const cc.std_string& path,
        const cc.bool is_sort = True):
    cdef list l = cc.qls(path, is_sort)
    return [ str(v) for v in l ]

def qls_all(const cc.std_string& path,
            const cc.bool is_folder_before_files = False,
            const cc.bool is_sort = True):
    cdef list l = cc.qls_all(path, is_folder_before_files, is_sort)
    return [ str(v) for v in l ]

def compute_crc32(const cc.std_string& path):
    return cc.compute_crc32(path)

def qload_datatable(const cc.std_string& path, const cc.bool is_par = False):
    return cc.qload_datatable(path, is_par)

def check_all_files_crc32_info(const cc.std_string& path):
    return cc.check_all_files_crc32_info(path)

