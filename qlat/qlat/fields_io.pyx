# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_base cimport FieldBase, SelectedFieldBase
from .field_selection cimport FieldSelection
from .field_types cimport FieldChar
from .selected_field_types cimport SelectedFieldChar

import qlat_utils as q

cache_fields_io = q.mk_cache("fields_io")

cdef class ShuffledFieldsWriter:

    # self.cdata

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __init__(self, str path not None, Coordinate new_size_node not None, cc.Bool is_append=False):
        self.xx.init(path, new_size_node.xx, is_append)

    def __del__(self):
        self.close()

    def close(self):
        self.xx.close()
        cache_fields_io.pop(id(self), None)

    def path(self):
        return self.xx.path

    def new_size_node(self):
        cdef Coordinate x = Coordinate()
        x.xx = self.xx.new_size_node
        return x

    def get_cache_sbs(self, FieldSelection fsel):
        if fsel is None:
            cache_fields_io.pop(id(self), None)
            return
        if id(self) in cache_fields_io:
            c_fsel, c_sbs = cache_fields_io[id(self)]
            if fsel is c_fsel:
                return c_sbs
        sbs = ShuffledBitSet(fsel, self.new_size_node())
        assert sbs.xx.fsel.n_elems == fsel.xx.n_elems
        assert cc.is_matching_fsel(sbs.xx.fsel, fsel.xx)
        cache_fields_io[id(self)] = (fsel, sbs,)
        return sbs

    def list(self):
        return cc.list_fields(self.xx)

    def has(self, const cc.std_string& fn):
        return cc.does_file_exist_sync_node(self.xx, fn)

    def __contains__(self, str fn):
        return self.has(fn)

    def flush(self):
        return cc.flush(self.xx)

    def write(self, str fn, obj):
        if isinstance(obj, (FieldBase, SelectedFieldBase,)):
            return obj.write_sfw_direct(self, fn)
        else:
            raise Exception("ShuffledFieldsWriter.write")

## --------------

cdef class ShuffledFieldsReader:

    # self.cdata

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __init__(self, str path not None, Coordinate new_size_node=None):
        if new_size_node is None:
            new_size_node = Coordinate()
        self.xx.init(path, new_size_node.xx)

    def __del__(self):
        self.close()

    def close(self):
        self.xx.close()
        cache_fields_io.pop(id(self), None)

    def path(self):
        return self.xx.path

    def new_size_node(self):
        cdef Coordinate x = Coordinate()
        x.xx = self.xx.new_size_node
        return x

    def get_cache_sbs(self, FieldSelection fsel):
        if fsel is None:
            cache_fields_io.pop(id(self), None)
            return
        if id(self) in cache_fields_io:
            c_fsel, c_sbs = cache_fields_io[id(self)]
            if fsel is c_fsel:
                return c_sbs
        sbs = ShuffledBitSet(fsel, self.new_size_node())
        assert sbs.xx.fsel.n_elems == fsel.xx.n_elems
        assert cc.is_matching_fsel(sbs.xx.fsel, fsel.xx)
        cache_fields_io[id(self)] = (fsel, sbs,)
        return sbs

    def list(self):
        return cc.list_fields(self.xx)

    def has_duplicates(self):
        return cc.has_duplicates(self.xx)

    def has(self, const cc.std_string& fn):
        return cc.does_file_exist_sync_node(self.xx, fn)

    def __contains__(self, str fn):
        return self.has(fn)

    def is_sparse_field(self, cc.std_string& fn):
        """
        return True if `fn in self` and is sparse field.
        """
        return cc.is_sparse_field_sync_node(self.xx, fn)

    def read_as_char(self, cc.std_string& fn):
        """
        return SelectedFieldChar or FieldChar.
        return None if `fn not in self`.
        """
        has_field = self.has(fn)
        if not has_field:
            return None
        is_sparse = self.is_sparse_field(fn)
        if is_sparse:
            obj = SelectedFieldChar(None)
        else:
            obj = FieldChar()
        total_bytes = obj.read_sfr_direct(self, fn)
        assert total_bytes > 0
        return obj

    def read(self, str fn, obj):
        """
        Can also read SelectedField obj with obj.fsel is None
        After reading, obj.fsel will be properly loaded.
        """
        if isinstance(obj, (FieldBase, SelectedFieldBase,)):
            return obj.read_sfr_direct(self, fn)
        else:
            raise Exception("ShuffledFieldsReader.read")

## --------------

cdef class ShuffledBitSet:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __init__(self, FieldSelection fsel not None, Coordinate new_size_node not None):
        self.xx = cc.mk_shuffled_bitset(fsel.xx, new_size_node.xx)

## --------------

@q.timer
def open_fields(str path, str mode, Coordinate new_size_node=None):
    """
    path can be the folder path or the 'geon-info.txt' path
    """
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if mode == "r":
        return ShuffledFieldsReader(path, new_size_node)
    elif mode == "w":
        assert new_size_node is not None
        return ShuffledFieldsWriter(path, new_size_node)
    elif mode == "a":
        if new_size_node is None:
            new_size_node = Coordinate()
        return ShuffledFieldsWriter(path, new_size_node, True)
    else:
        raise Exception("open_fields")

@q.timer
def list_fields(str path, Coordinate new_size_node=None):
    cdef ShuffledFieldsReader sfr = open_fields(path, "r", new_size_node)
    cdef list fns = sfr.list()
    sfr.close()
    return fns

@q.timer
def fields_build_index(str path, Coordinate new_size_node=None):
    cdef ShuffledFieldsReader sfr = open_fields(path, "r", new_size_node)
    cc.fields_build_index(sfr.xx)
    sfr.close()

@q.timer
def fields_has_duplicates(str path, Coordinate new_size_node=None):
    cdef ShuffledFieldsReader sfr = open_fields(path, "r", new_size_node)
    ret = sfr.has_duplicates();
    sfr.close()
    return ret

@q.timer
def properly_truncate_fields(str path, cc.Bool is_check_all=False, cc.Bool is_only_check=False, Coordinate new_size_node=None):
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if new_size_node is None:
        new_size_node = Coordinate()
    return cc.properly_truncate_fields_sync_node(path, is_check_all, is_only_check, new_size_node.xx)

@q.timer
def truncate_fields(str path, list fns_keep, Coordinate new_size_node=None):
    """
    fns_keep is the list of fields that need to keep
    fns_keep needs to be in the same order as the data is stored in path
    all fns_keep must be already in the path
    fns_keep can be empty list
    """
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if new_size_node is None:
        new_size_node = Coordinate()
    ret = cc.truncate_fields_sync_node(path, fns_keep, new_size_node.xx)
    if ret != 0:
        raise Exception(f"truncate_fields: error {ret}")

@q.timer
def check_fields(str path, cc.Bool is_check_all=True, Coordinate new_size_node=None):
    """
    return list of field that is stored successful
    """
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    is_only_check = True
    return properly_truncate_fields(path, is_check_all, is_only_check, new_size_node)

@q.timer
def check_compressed_eigen_vectors(str path):
    """
    return bool value suggest whether the data can be read successfully
    return True is the data has problem
    return False if the data is ok
    """
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    return cc.check_compressed_eigen_vectors(path)

@q.timer
def eigen_system_repartition(Coordinate new_size_node, str path, str path_new=""):
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if path_new[-14:] == "/geon-info.txt":
        path_new = path_new[:-14]
    return cc.eigen_system_repartition(new_size_node.xx, path, path_new)

## --------------

def show_all_shuffled_fields_writer():
    return cc.show_all_shuffled_fields_writer()
