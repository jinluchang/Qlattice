# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc

from .timer import timer

### ----------------------------------------------------------

cdef class QFile:

    def __init__(self, *,
            str path=None, str mode=None,
            QFile qfile=None, cc.Long q_offset_start=0, cc.Long q_offset_end=-1,
            ):
        """
        QFile(path=path, mode=mode)
        QFile(qfile=path, q_offset_start=q_offset_start, q_offset_end=q_offset_end)
        """
        if path is not None and qfile is None:
            self.xx.init(path, cc.read_qfile_mode(mode))
        elif path is None and qfile is not None:
            self.xx.init(qfile.xx, q_offset_start, q_offset_end)
        else:
            assert path is None
            assert qfile is None

    def __imatmul__(self, QFile v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data=True):
        cdef QFile x = type(self)()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def path(self):
        return self.xx.path()

    def mode(self):
        return cc.show(self.xx.mode())

    def close(self):
        return self.xx.close()

    def null(self):
        return self.xx.null()

    def eof(self):
        return cc.qfeof(self.xx)

    def tell(self):
        return cc.qftell(self.xx)

    def flush(self):
        return cc.qfflush(self.xx)

    def seek_set(self, cc.Long q_offset):
        return cc.qfseek_set(self.xx, q_offset)

    def seek_end(self, cc.Long q_offset):
        return cc.qfseek_end(self.xx, q_offset)

    def seek_cur(self, cc.Long q_offset):
        return cc.qfseek_cur(self.xx, q_offset)

    def size(self):
        return cc.qfile_size(self.xx)

    def remaining_size(self):
        return cc.qfile_remaining_size(self.xx)

    def getline(self):
        return cc.qgetline(self.xx)

    def getlines(self):
        return cc.qgetlines(self.xx)

    def getlines(self):
        return cc.qgetlines(self.xx)

    def qcat(self):
        return <str>cc.qcat(self.xx)

    def qcat_bytes(self):
        return <bytes>cc.qcat(self.xx)

    def write(self, object content):
        cdef QFile qfile
        if isinstance(content, QFile):
            qfile = content
            return cc.write_from_qfile(self.xx, qfile.xx)
        elif isinstance(content, (str, bytes,)):
            return cc.qwrite_data(<cc.std_string>content, self.xx)
        elif isinstance(content, (list, tuple,)):
            return cc.qwrite_data(<cc.std_vector[cc.std_string]>content, self.xx)
        else:
            raise Exception(f"write: {type(content)}")

    def compute_crc32(self):
        return cc.compute_crc32(self.xx)

### ----------------------------------------------------------

cdef class QarFile:

    def __init__(self, *, str path=None, str mode=None):
        """
        QarFile(path=path, mode=mode)
        """
        if path is not None:
            self.xx.init(path, cc.read_qfile_mode(mode))
        else:
            assert path is None

    def __imatmul__(self, QarFile v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data=True):
        cdef QarFile x = type(self)()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def path(self):
        return self.xx.path

    def mode(self):
        return cc.show(self.xx.mode)

    def close(self):
        return self.xx.close()

    def null(self):
        return self.xx.null()

    def list(self):
        return cc.list(self.xx)

    def has_regular_file(self, const cc.std_string& fn):
        return cc.has_regular_file(self.xx, fn)

    def has(self, const cc.std_string& fn):
        return cc.has(self.xx, fn)

    def read(self, const cc.std_string& fn):
        """
        return qfile
        """
        cdef QFile qfile = QFile()
        qfile.xx = cc.read(self.xx, fn)
        return qfile

    def read_data(self, const cc.std_string& fn):
        return <str>cc.read_data(self.xx, fn)

    def read_data_bytes(self, const cc.std_string& fn):
        return <bytes>cc.read_data(self.xx, fn)

    def read_info(self, const cc.std_string& fn):
        return cc.read_info(self.xx, fn)

    def verify_index(self):
        return cc.verify_index(self.xx)

    def write(self, const cc.std_string& fn, const cc.std_string& info, object data):
        cdef QFile qfile
        if isinstance(data, QFile):
            qfile = data
            return cc.write_from_qfile(self.xx, fn, info, qfile.xx)
        elif isinstance(data, (str, bytes,)):
            return cc.write_from_data(self.xx, fn, info, <cc.std_string>data)
        elif isinstance(data, (list, tuple,)):
            return cc.write_from_data(self.xx, fn, info, <cc.std_vector[cc.std_string]>data)
        else:
            raise Exception(f"write: {type(data)}")

    def show_index(self):
        return cc.show_qar_index(self.xx)

    def save_index(self, const cc.std_string& fn):
        return cc.save_qar_index(self.xx, fn)

    def parse_index(self, const cc.std_string& qar_index_content):
        return cc.parse_qar_index(self.xx, qar_index_content)

    def load_index(self, const cc.std_string& fn):
        return cc.load_qar_index(self.xx, fn)

### ----------------------------------------------------------

def get_qar_multi_vol_max_size():
    """
    Parameter controls the size of a single `qar` file in number of bytes. Note, `qar` never splits a single file into multiple `qar` volume.
    """
    return cc.get_qar_multi_vol_max_size()

def set_qar_multi_vol_max_size(size=None):
    """
    Parameter controls the size of a single `qar` file in number of bytes. Note, `qar` never splits a single file into multiple `qar` volume.
    """
    if size is None:
        size = cc.get_qar_multi_vol_max_size_default()
    assert isinstance(size, int)
    cdef cc.Long* p_size = &cc.get_qar_multi_vol_max_size()
    p_size[0] = size
    assert cc.get_qar_multi_vol_max_size() == size

def clean_up_qfile_map():
    return cc.clean_up_qfile_map()

def show_all_qfile():
    return cc.show_all_qfile()

### ----------------------------------------------------------

def properly_truncate_qar_file(const cc.std_string& path):
    return cc.properly_truncate_qar_file(path)

### ----------------------------------------------------------

def does_regular_file_exist_qar(const cc.std_string& path):
    return cc.does_regular_file_exist_qar(path)

def does_file_exist_qar(const cc.std_string& path):
    return cc.does_file_exist_qar(path)

@timer
def qar_build_index(const cc.std_string& path_qar):
    """
    create "path_qar.idx" file
    """
    cc.qar_build_index(path_qar)

@timer
def qar_create(const cc.std_string& path_qar, const cc.std_string& path_folder,
               *, const cc.bool is_remove_folder_after=False):
    return cc.qar_create(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract(const cc.std_string& path_qar, const cc.std_string& path_folder,
                *, const cc.bool is_remove_qar_after=False):
    return cc.qar_extract(path_qar, path_folder, is_remove_qar_after)

@timer
def qcopy_file(const cc.std_string& path_src, const cc.std_string& path_dst):
    return cc.qcopy_file(path_src, path_dst)

@timer
def list_qar(const cc.std_string& path_qar):
    cdef list l = cc.list_qar(path_qar)
    return [ <str>fn for fn in l ]

### ----------------------------------------------------------

def qcat(const cc.std_string& path):
    """Return contents of file as `str`"""
    return <str>cc.qcat(path)

def qcat_bytes(const cc.std_string& path):
    """Return contents of file as `bytes`"""
    return <bytes>cc.qcat(path)

def qtouch(const cc.std_string& path, object content=None):
    if content is None:
        return cc.qtouch(path)
    elif isinstance(content, (str, bytes,)):
        return cc.qtouch(path, <cc.std_string>content)
    elif isinstance(content, (list, tuple,)):
        return cc.qtouch(path, <cc.std_vector[cc.std_string]>content)
    else:
        raise Exception(f"qtouch: {type(content)}")

def qappend(const cc.std_string& path, object content):
    if isinstance(content, (str, bytes,)):
        return cc.qappend(path, <cc.std_string>content)
    elif isinstance(content, (list, tuple,)):
        return cc.qappend(path, <cc.std_vector[cc.std_string]>content)
    else:
        raise Exception(f"qappend: {type(content)}")

def qload_datatable(const cc.std_string& path, const cc.bool is_par=False):
    return cc.qload_datatable(path, is_par)

def compute_crc32(const cc.std_string& path):
    return cc.compute_crc32(path)

### ----------------------------------------------------------

@timer
def qar_build_index_info(const cc.std_string& path_qar):
    """
    create "path_qar.idx" file (only on node 0)
    """
    cc.qar_build_index_info(path_qar)

@timer
def qar_create_info(const cc.std_string& path_qar, const cc.std_string& path_folder,
               *, const cc.bool is_remove_folder_after=False):
    return cc.qar_create_info(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract_info(const cc.std_string& path_qar, const cc.std_string& path_folder,
                *, const cc.bool is_remove_qar_after=False):
    return cc.qar_extract_info(path_qar, path_folder, is_remove_qar_after)

@timer
def qcopy_file_info(const cc.std_string& path_src, const cc.std_string& path_dst):
    return cc.qcopy_file_info(path_src, path_dst)

def qtouch_info(const cc.std_string& path, object content=None):
    if content is None:
        return cc.qtouch_info(path)
    elif isinstance(content, (str, bytes,)):
        return cc.qtouch_info(path, <cc.std_string>content)
    elif isinstance(content, (list, tuple,)):
        return cc.qtouch_info(path, <cc.std_vector[cc.std_string]>content)
    else:
        raise Exception(f"qtouch: {type(content)}")

def qappend_info(const cc.std_string& path, object content):
    if isinstance(content, (str, bytes,)):
        return cc.qappend_info(path, <cc.std_string>content)
    elif isinstance(content, (list, tuple,)):
        return cc.qappend_info(path, <cc.std_vector[cc.std_string]>content)
    else:
        raise Exception(f"qappend: {type(content)}")

def check_all_files_crc32_info(const cc.std_string& path):
    return cc.check_all_files_crc32_info(path)

### ----------------------------------------------------------

def does_regular_file_exist_qar_sync_node(const cc.std_string& path):
    return cc.does_regular_file_exist_qar_sync_node(path)

def does_file_exist_qar_sync_node(const cc.std_string& path):
    return cc.does_file_exist_qar_sync_node(path)

@timer
def qar_create_sync_node(const cc.std_string& path_qar, const cc.std_string& path_folder,
        *, const cc.bool is_remove_folder_after=False):
    return cc.qar_create_sync_node(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract_sync_node(const cc.std_string& path_qar, const cc.std_string& path_folder,
                *, const cc.bool is_remove_qar_after=False):
    return cc.qar_extract_sync_node(path_qar, path_folder, is_remove_qar_after)

@timer
def qcopy_file_sync_node(const cc.std_string& path_src, const cc.std_string& path_dst):
    return cc.qcopy_file_sync_node(path_src, path_dst)

def qcat_sync_node(const cc.std_string& path):
    """Return contents of file as `str`"""
    return <str>cc.qcat_sync_node(path)

def qcat_bytes_sync_node(const cc.std_string& path):
    """Return contents of file as `bytes`"""
    return <bytes>cc.qcat_sync_node(path)

def qload_datatable_sync_node(const cc.std_string& path, const cc.bool is_par=False):
    return cc.qload_datatable_sync_node(path, is_par)

### ----------------------------------------------------------
