# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc

from .timer import timer, get_id_node

### ----------------------------------------------------------

class Gobble:
    def __getattr__(self, item):
        return self
    def __call__(self, *args, **kwargs):
        return self

cdef class QFile:

    def __init__(self, *,
            str ftype=None, str path=None, str mode=None, object content=None,
            QFile qfile=None, cc.Long q_offset_start=0, cc.Long q_offset_end=-1,
            ):
        """
        Always use kwargs.
        #
        `QFile(path, mode)`
        `QFile(ftype, path, mode)`
        `QFile(ftype, path, mode, content)`
        `QFile(qfile, q_offset_start, q_offset_end)`
        #
        `mode in [ "r", "w", "a", ]`
        """
        cdef cc.QFileType ftype_v
        cdef cc.std_string path_v
        cdef cc.QFileMode mode_v
        cdef cc.std_string content_v
        cdef str ftype_default
        if ftype is not None:
            ftype_v = cc.read_qfile_type(ftype)
        else:
            ftype_default = "CFile"
            ftype_v = cc.read_qfile_type(ftype_default)
        if path is not None:
            path_v = path
        if mode is not None:
            mode_v = cc.read_qfile_mode(mode)
        else:
            mode_default = "r"
            mode_v = cc.read_qfile_mode(mode_default)
        if content is not None:
            content_v = content
        if content is not None:
            assert ftype is not None
            assert path is not None
            assert mode is not None
            assert qfile is None
            self.xx.init(ftype_v, path_v, mode_v, content_v)
        elif ftype is not None:
            assert path is not None
            assert mode is not None
            assert content is None
            assert qfile is None
            self.xx.init(ftype_v, path_v, mode_v)
        elif path is not None:
            assert ftype is None
            assert mode is not None
            assert content is None
            assert qfile is None
            self.xx.init(ftype_v, path_v, mode_v)
        elif qfile is not None:
            assert ftype is None
            assert path is None
            assert mode is None
            assert content is None
            self.xx.init(qfile.xx, q_offset_start, q_offset_end)
        else:
            assert ftype is None
            assert path is None
            assert mode is None
            assert content is None
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

    def content(self):
        cdef cc.std_string x = self.xx.content()
        return <str>x;

    def content_bytes(self):
        cdef cc.std_string x = self.xx.content()
        return <bytes>x;

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

    def write(self, object data):
        cdef QFile qfile
        if isinstance(data, QFile):
            qfile = data
            return cc.write_from_qfile(self.xx, qfile.xx)
        elif isinstance(data, (str, bytes,)):
            return cc.qwrite_data(<cc.std_string>data, self.xx)
        elif isinstance(data, (list, tuple,)):
            return cc.qwrite_data(<cc.std_vector[cc.std_string]>data, self.xx)
        else:
            raise Exception(f"write: {type(data)}")

    def compute_crc32(self):
        return cc.compute_crc32(self.xx)

### ----------------------------------------------------------

def open_qfile(const cc.std_string& path, const cc.std_string& mode):
    """
    Call `cc.qfopen`
    """
    cdef cc.std_string path_v = path
    cdef cc.QFileMode mode_v = cc.read_qfile_mode(mode)
    cdef QFile qfile = QFile()
    qfile.xx = cc.qfopen(path_v, mode_v);
    return qfile

def open_qfile_str(const cc.std_string& path, const cc.std_string& mode, object content=None):
    """
    Call `cc.qfopen`
    """
    cdef str ftype = "String"
    cdef cc.QFileType ftype_v = cc.read_qfile_type(ftype)
    cdef cc.std_string path_v = path
    cdef cc.QFileMode mode_v = cc.read_qfile_mode(mode)
    cdef cc.std_string content_v
    cdef QFile qfile = QFile()
    if content is None:
        qfile.xx = cc.qfopen(ftype_v, path_v, mode_v);
    else:
        content_v = content
        qfile.xx = cc.qfopen(ftype_v, path_v, mode_v, content_v);
    return qfile

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
            assert mode is None

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

    def flush(self):
        return self.xx.flush()

    def list(self):
        return cc.list(self.xx)

    def has_regular_file(self, const cc.std_string& fn):
        return cc.has_regular_file(self.xx, fn)

    def has(self, const cc.std_string& fn):
        return cc.has(self.xx, fn)

    def __contains__(self, str fn):
        return self.has(fn)

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

    def write(self, const cc.std_string& fn, const cc.std_string& info, object data, *, cc.Bool skip_if_exist=False):
        cdef QFile qfile
        if skip_if_exist:
            if self.has(fn):
                return 0
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

    def read_index(self, const cc.std_string& qar_index_content):
        return cc.read_qar_index(self.xx, qar_index_content)

    def index_size_saved(self):
        cdef cc.Long size = self.xx.qar_index_size_saved
        return size

    def index_size(self):
        cdef cc.Long size = self.xx.index_size()
        return size

    def save_index(self, cc.Long max_diff=0):
        self.xx.save_index(max_diff)

### ----------------------------------------------------------

@timer
def open_qar(const cc.std_string& path, const cc.std_string& mode):
    """
    Call QarFile(path, mode) with kwargs
    """
    cdef QarFile qar = QarFile(path=path, mode=mode)
    return qar

def open_qar_info(*args, **kwargs):
    """
    Call `open_qar` with same arguments if q.get_id_node() == 0.
    Otherwise return Gobble(), which does nothing for any method and returns it self.
    """
    if get_id_node() == 0:
        return open_qar(*args, **kwargs)
    else:
        return Gobble()

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

def qload_datatable(const cc.std_string& path, const cc.Bool is_par=False):
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
               *, const cc.Bool is_remove_folder_after=False):
    return cc.qar_create_info(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract_info(const cc.std_string& path_qar, const cc.std_string& path_folder,
                *, const cc.Bool is_remove_qar_after=False):
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
        *, const cc.Bool is_remove_folder_after=False):
    return cc.qar_create_sync_node(path_qar, path_folder, is_remove_folder_after)

@timer
def qar_extract_sync_node(const cc.std_string& path_qar, const cc.std_string& path_folder,
                *, const cc.Bool is_remove_qar_after=False):
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

def qload_datatable_sync_node(const cc.std_string& path, const cc.Bool is_par=False):
    return cc.qload_datatable_sync_node(path, is_par)

### ----------------------------------------------------------
