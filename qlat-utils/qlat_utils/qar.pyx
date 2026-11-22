# cython: binding=True, embedsignature=True

"""
Module ``qlat_utils.qar``
=========================\n
QAR archive format, QFile abstraction, and file I/O utilities.\n
Documentation: ``docs/qlat-utils/qlat_qar.md``\n
.. note:: Update the documentation when updating this source file.
"""

from . cimport everything as cc

class q:
    from .timer import (
        timer,
        get_id_node,
    )

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
            ftype_v = cc.read_qfile_type(cc.as_bytes(ftype))
        else:
            ftype_default = "CFile"
            ftype_v = cc.read_qfile_type(cc.as_bytes(ftype_default))
        if path is not None:
            path_v = cc.as_bytes(path)
        if mode is not None:
            mode_v = cc.read_qfile_mode(cc.as_bytes(mode))
        else:
            mode_default = "r"
            mode_v = cc.read_qfile_mode(cc.as_bytes(mode_default))
        if content is not None:
            content_v = cc.as_bytes(content)
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
        return cc.as_str_std(self.xx.path())

    def mode(self):
        return cc.as_str_std(cc.show(self.xx.mode()))

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
        return cc.as_str_std(x)

    def content_bytes(self):
        cdef cc.std_string x = self.xx.content()
        return <bytes>x;

    def size(self):
        return cc.qfile_size(self.xx)

    def remaining_size(self):
        return cc.qfile_remaining_size(self.xx)

    def getline(self):
        return cc.as_str_std(cc.qgetline(self.xx))

    def getlines(self):
        return [cc.as_str_std(v) for v in cc.qgetlines(self.xx)]

    def qcat(self):
        return cc.as_str_std(cc.qcat(self.xx))

    def qcat_bytes(self):
        return <bytes>cc.qcat(self.xx)

    def write(self, object data):
        cdef QFile qfile
        cdef cc.std_vector[cc.std_string] data_v
        if isinstance(data, QFile):
            qfile = data
            return cc.write_from_qfile(self.xx, qfile.xx)
        elif isinstance(data, (str, bytes,)):
            return cc.qwrite_data(cc.as_bytes(data), self.xx)
        elif isinstance(data, (list, tuple,)):
            data_v = [cc.as_bytes(x) for x in data]
            return cc.qwrite_data(data_v, self.xx)
        else:
            raise Exception(f"write: {type(data)}")

    def compute_crc32(self):
        return cc.compute_crc32(self.xx)

### ----------------------------------------------------------

def open_qfile(path, mode):
    """
    Call `cc.qfopen`
    """
    cdef bytes path_b = cc.as_bytes(path)
    cdef bytes mode_b = cc.as_bytes(mode)
    cdef cc.std_string path_v = path_b
    cdef cc.QFileMode mode_v = cc.read_qfile_mode(mode_b)
    cdef QFile qfile = QFile()
    qfile.xx = cc.qfopen(path_v, mode_v);
    return qfile

def open_qfile_str(path, mode, object content=None):
    """
    Call `cc.qfopen`
    """
    cdef bytes path_b = cc.as_bytes(path)
    cdef bytes mode_b = cc.as_bytes(mode)
    cdef bytes ftype = b"String"
    cdef cc.QFileType ftype_v = cc.read_qfile_type(ftype)
    cdef cc.std_string path_v = path_b
    cdef cc.QFileMode mode_v = cc.read_qfile_mode(mode_b)
    cdef cc.std_string content_v
    cdef QFile qfile = QFile()
    if content is None:
        qfile.xx = cc.qfopen(ftype_v, path_v, mode_v);
    else:
        content_v = cc.as_bytes(content)
        qfile.xx = cc.qfopen(ftype_v, path_v, mode_v, content_v);
    return qfile

### ----------------------------------------------------------

cdef class QarFile:

    def __init__(self, *, str path=None, str mode=None):
        """
        QarFile(path=path, mode=mode)\n
        Parameters are keyword-only.  Use::\n
            QarFile(path="/path/to/file.qar", mode="w")\n
        or the ``open_qar()`` convenience wrapper which accepts positional args.
        """
        if path is not None:
            self.xx.init(cc.as_bytes(path), cc.read_qfile_mode(cc.as_bytes(mode)))
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
        return cc.as_str_std(self.xx.path)

    def mode(self):
        return cc.as_str_std(cc.show(self.xx.mode))

    def close(self):
        return self.xx.close()

    def null(self):
        return self.xx.null()

    def flush(self):
        return self.xx.flush()

    def list(self):
        """
        List entries in this QAR archive.
        #
        When the QAR was opened via ``open_qar_info``, this method returns a
        ``Gobble`` on non-rank-0 nodes.  Broadcast the result from rank 0 with
        ``q.bcast_py`` to synchronize the listing across MPI ranks.
        """
        return [cc.as_str_std(v) for v in cc.list(self.xx)]

    def has_regular_file(self, fn):
        """
        Check if a regular file entry exists in this QAR archive.
        #
        When the QAR was opened via ``open_qar_info``, this method returns a
        truthy ``Gobble`` on non-rank-0 nodes (because ``open_qar_info`` returns
        ``Gobble()`` on rank > 0).  To avoid a control-flow desync between MPI
        ranks, broadcast the result from rank 0 with ``q.bcast_py``::
            #
            has_file = q.bcast_py(qar_ws.has_regular_file(fn))
        """
        cdef bytes fn_b = cc.as_bytes(fn)
        return cc.has_regular_file(self.xx, fn_b)

    def has(self, fn):
        """
        Check if an entry exists in this QAR archive (regular file or otherwise).
        #
        When the QAR was opened via ``open_qar_info``, this method returns a
        truthy ``Gobble`` on non-rank-0 nodes.  Broadcast the result from rank 0
        with ``q.bcast_py`` to avoid a control-flow desync between MPI ranks::
            #
            has_entry = q.bcast_py(qar_ws.has(fn))
        """
        cdef bytes fn_b = cc.as_bytes(fn)
        return cc.has(self.xx, fn_b)

    def __contains__(self, str fn):
        """
        ``fn in qar_ws`` — delegates to :meth:`has`.  The same ``Gobble`` / ``bcast_py``
        caveat applies when the QAR was opened via ``open_qar_info``.
        """
        return self.has(fn)

    def read(self, fn):
        """
        return qfile
        """
        cdef QFile qfile = QFile()
        cdef bytes fn_b = cc.as_bytes(fn)
        qfile.xx = cc.read(self.xx, fn_b)
        return qfile

    def read_data(self, fn):
        """
        Read entry data as a Python ``str`` (UTF-8 decoded).\n
        Only suitable for text content.  For binary data use
        :meth:`read_data_bytes` instead — this method will raise
        ``UnicodeDecodeError`` on non-UTF-8 bytes.
        """
        cdef bytes fn_b = cc.as_bytes(fn)
        return cc.as_str_std(cc.read_data(self.xx, fn_b))

    def read_data_bytes(self, fn):
        """
        Read entry data as raw ``bytes``.\n
        Safe for both text and binary content.  Prefer this over
        :meth:`read_data` unless you know the entry is plain text.
        """
        cdef bytes fn_b = cc.as_bytes(fn)
        return <bytes>cc.read_data(self.xx, fn_b)

    def read_info(self, fn):
        cdef bytes fn_b = cc.as_bytes(fn)
        return cc.read_info(self.xx, fn_b)

    def verify_index(self):
        return cc.verify_index(self.xx)

    def write(self, fn, info, object data, *, cc.Bool skip_if_exist=False):
        cdef bytes fn_b = cc.as_bytes(fn)
        cdef bytes info_b = cc.as_bytes(info)
        cdef cc.std_vector[cc.std_string] data_v
        cdef QFile qfile
        if skip_if_exist:
            if self.has(fn):
                return 0
        if isinstance(data, QFile):
            qfile = data
            return cc.write_from_qfile(self.xx, fn_b, info_b, qfile.xx)
        elif isinstance(data, (str, bytes,)):
            return cc.write_from_data(self.xx, fn_b, info_b, cc.as_bytes(data))
        elif isinstance(data, (list, tuple,)):
            data_v = [cc.as_bytes(x) for x in data]
            return cc.write_from_data(self.xx, fn_b, info_b, data_v)
        else:
            raise Exception(f"write: {type(data)}")

    def show_index(self):
        return cc.show_qar_index(self.xx)

    def read_index(self, qar_index_content):
        cdef bytes content_b = cc.as_bytes(qar_index_content)
        return cc.read_qar_index(self.xx, content_b)

    def index_size_saved(self):
        cdef cc.Long size = self.xx.qar_index_size_saved
        return size

    def index_size(self):
        cdef cc.Long size = self.xx.index_size()
        return size

    def save_index(self, cc.Long max_diff=0):
        self.xx.save_index(max_diff)

### ----------------------------------------------------------

@q.timer
def open_qar(path, mode):
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
    if q.get_id_node() == 0:
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
    return [cc.as_str_std(v) for v in cc.show_all_qfile()]

### ----------------------------------------------------------

def properly_truncate_qar_file(path):
    cdef bytes path_b = cc.as_bytes(path)
    return [cc.as_str_std(v) for v in cc.properly_truncate_qar_file(path_b)]

### ----------------------------------------------------------

def does_regular_file_exist_qar(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.does_regular_file_exist_qar(path_b)

def does_file_exist_qar(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.does_file_exist_qar(path_b)

@q.timer
def qar_build_index(path_qar):
    """
    create "path_qar.idx" file
    """
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    cc.qar_build_index(path_qar_b)

@q.timer
def qar_create(path_qar, path_folder,
               *, const cc.bool is_remove_folder_after=False):
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    cdef bytes path_folder_b = cc.as_bytes(path_folder)
    return cc.qar_create(path_qar_b, path_folder_b, is_remove_folder_after)

@q.timer
def qar_extract(path_qar, path_folder,
                *, const cc.bool is_remove_qar_after=False):
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    cdef bytes path_folder_b = cc.as_bytes(path_folder)
    return cc.qar_extract(path_qar_b, path_folder_b, is_remove_qar_after)

@q.timer
def qcopy_file(path_src, path_dst):
    cdef bytes path_src_b = cc.as_bytes(path_src)
    cdef bytes path_dst_b = cc.as_bytes(path_dst)
    return cc.qcopy_file(path_src_b, path_dst_b)

@q.timer
def list_qar(path_qar):
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    return [cc.as_str_std(v) for v in cc.list_qar(path_qar_b)]

### ----------------------------------------------------------

def qcat(path):
    """Return contents of file as `str`"""
    cdef bytes path_b = cc.as_bytes(path)
    return cc.as_str_std(cc.qcat(path_b))

def qcat_bytes(path):
    """Return contents of file as `bytes`"""
    cdef bytes path_b = cc.as_bytes(path)
    return <bytes>cc.qcat(path_b)

def qtouch(path, object content=None):
    cdef bytes path_b = cc.as_bytes(path)
    cdef cc.std_vector[cc.std_string] content_v
    if content is None:
        return cc.qtouch(path_b)
    elif isinstance(content, (str, bytes,)):
        return cc.qtouch(path_b, cc.as_bytes(content))
    elif isinstance(content, (list, tuple,)):
        content_v = [cc.as_bytes(x) for x in content]
        return cc.qtouch(path_b, content_v)
    else:
        raise Exception(f"qtouch: {type(content)}")

def qappend(path, object content):
    cdef bytes path_b = cc.as_bytes(path)
    cdef cc.std_vector[cc.std_string] content_v
    if isinstance(content, (str, bytes,)):
        return cc.qappend(path_b, cc.as_bytes(content))
    elif isinstance(content, (list, tuple,)):
        content_v = [cc.as_bytes(x) for x in content]
        return cc.qappend(path_b, content_v)
    else:
        raise Exception(f"qappend: {type(content)}")

def qload_datatable(path, const cc.Bool is_par=False):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qload_datatable(path_b, is_par)

def compute_crc32(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.compute_crc32(path_b)

### ----------------------------------------------------------

@q.timer
def qar_build_index_info(path_qar):
    """
    create "path_qar.idx" file (only on node 0)
    """
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    cc.qar_build_index_info(path_qar_b)

@q.timer
def qar_create_info(path_qar, path_folder,
               *, const cc.Bool is_remove_folder_after=False):
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    cdef bytes path_folder_b = cc.as_bytes(path_folder)
    return cc.qar_create_info(path_qar_b, path_folder_b, is_remove_folder_after)

@q.timer
def qar_extract_info(path_qar, path_folder,
                *, const cc.Bool is_remove_qar_after=False):
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    cdef bytes path_folder_b = cc.as_bytes(path_folder)
    return cc.qar_extract_info(path_qar_b, path_folder_b, is_remove_qar_after)

@q.timer
def qcopy_file_info(path_src, path_dst):
    cdef bytes path_src_b = cc.as_bytes(path_src)
    cdef bytes path_dst_b = cc.as_bytes(path_dst)
    return cc.qcopy_file_info(path_src_b, path_dst_b)

def qtouch_info(path, object content=None):
    cdef bytes path_b = cc.as_bytes(path)
    cdef cc.std_vector[cc.std_string] content_v
    if content is None:
        return cc.qtouch_info(path_b)
    elif isinstance(content, (str, bytes,)):
        return cc.qtouch_info(path_b, cc.as_bytes(content))
    elif isinstance(content, (list, tuple,)):
        content_v = [cc.as_bytes(x) for x in content]
        return cc.qtouch_info(path_b, content_v)
    else:
        raise Exception(f"qtouch: {type(content)}")

def qappend_info(path, object content):
    cdef bytes path_b = cc.as_bytes(path)
    cdef cc.std_vector[cc.std_string] content_v
    if isinstance(content, (str, bytes,)):
        return cc.qappend_info(path_b, cc.as_bytes(content))
    elif isinstance(content, (list, tuple,)):
        content_v = [cc.as_bytes(x) for x in content]
        return cc.qappend_info(path_b, content_v)
    else:
        raise Exception(f"qappend: {type(content)}")

def check_all_files_crc32_info(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.check_all_files_crc32_info(path_b)

### ----------------------------------------------------------

def does_regular_file_exist_qar_sync_node(path):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.does_regular_file_exist_qar_sync_node(path_b)

def does_file_exist_qar_sync_node(path):
    """
    Return True if file or directory exists.
    """
    cdef bytes path_b = cc.as_bytes(path)
    return cc.does_file_exist_qar_sync_node(path_b)

@q.timer
def qar_create_sync_node(
    path_qar, path_folder,
    *, const cc.Bool is_remove_folder_after=False,
):
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    cdef bytes path_folder_b = cc.as_bytes(path_folder)
    return cc.qar_create_sync_node(path_qar_b, path_folder_b, is_remove_folder_after)

@q.timer
def qar_extract_sync_node(path_qar, path_folder,
                *, const cc.Bool is_remove_qar_after=False):
    cdef bytes path_qar_b = cc.as_bytes(path_qar)
    cdef bytes path_folder_b = cc.as_bytes(path_folder)
    return cc.qar_extract_sync_node(path_qar_b, path_folder_b, is_remove_qar_after)

@q.timer
def qcopy_file_sync_node(path_src, path_dst):
    cdef bytes path_src_b = cc.as_bytes(path_src)
    cdef bytes path_dst_b = cc.as_bytes(path_dst)
    return cc.qcopy_file_sync_node(path_src_b, path_dst_b)

def qcat_sync_node(path):
    """Return contents of file as `str`"""
    cdef bytes path_b = cc.as_bytes(path)
    return cc.as_str_std(cc.qcat_sync_node(path_b))

def qcat_bytes_sync_node(path):
    """Return contents of file as `bytes`"""
    cdef bytes path_b = cc.as_bytes(path)
    return <bytes>cc.qcat_sync_node(path_b)

def qload_datatable_sync_node(path, const cc.Bool is_par=False):
    cdef bytes path_b = cc.as_bytes(path)
    return cc.qload_datatable_sync_node(path_b, is_par)

### ----------------------------------------------------------
