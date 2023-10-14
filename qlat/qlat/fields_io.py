from qlat_utils import *
import qlat.c as c

from qlat.field import *
from qlat.field_selection import *
from qlat.selected_field import *
from qlat.utils_io import *

cache_fields_io = mk_cache("fields_io")

class ShuffledFieldsWriter:

    # self.cdata

    def __init__(self, path, new_size_node, is_append = False):
        assert isinstance(path, str)
        assert isinstance(is_append , bool)
        self.cdata = c.mk_sfw(path, new_size_node, is_append)

    def close(self):
        if self.cdata is not None:
            c.free_sfw(self)
        self.cdata = None
        cache_fields_io.pop(id(self), None)

    def __del__(self):
        self.close()

    def new_size_node(self):
        return c.get_new_size_node_sfw(self)

    def get_cache_sbs(self, fsel):
        if id(self) in cache_fields_io:
            [ c_fsel, c_sbs ] = cache_fields_io[id(self)]
            if fsel is c_fsel:
                return c_sbs
        sbs = ShuffledBitSet(fsel, self.new_size_node())
        cache_fields_io[id(self)] = [ fsel, sbs, ]
        return sbs

    def write(self, fn, obj):
        assert isinstance(fn, str)
        if isinstance(obj, FieldBase):
            return c.write_sfw_field(self, fn, obj)
        elif isinstance(obj, SelectedFieldBase):
            return c.write_sfw_sfield(self, fn, obj, self.get_cache_sbs(obj.fsel))
        else:
            raise Exception("ShuffledFieldsWriter.save")

    def flush(self):
        return c.flush_sfw(self)

class ShuffledFieldsReader:

    # self.cdata
    # self.tags

    def __init__(self, path, new_size_node = None):
        assert isinstance(path, str)
        if new_size_node is None:
            self.cdata = c.mk_sfr(path)
        else:
            self.cdata = c.mk_sfr(path, new_size_node)
        self.tags = None

    def close(self):
        if self.cdata is not None:
            c.free_sfr(self)
        self.cdata = None
        self.tags = None
        cache_fields_io.pop(id(self), None)

    def __del__(self):
        self.close()

    def new_size_node(self):
        return c.get_new_size_node_sfr(self)

    def get_cache_sbs(self, fsel):
        if id(self) in cache_fields_io:
            [ c_fsel, c_sbs ] = cache_fields_io[id(self)]
            if fsel is c_fsel:
                return c_sbs
        sbs = ShuffledBitSet(fsel, self.new_size_node())
        cache_fields_io[id(self)] = [ fsel, sbs, ]
        return sbs

    def read(self, fn, obj):
        # Can also read SelectedField obj with obj.fsel is None
        # After reading, obj.fsel will be properly loaded.
        assert isinstance(fn, str)
        if isinstance(obj, FieldBase):
            return c.read_sfr_field(self, fn, obj)
        elif isinstance(obj, SelectedFieldBase):
            fsel = obj.fsel
            if fsel is None:
                if obj.view_count > 0:
                    raise ValueError("can't re-init while being viewed")
                obj.fsel = FieldSelection()
                return c.read_sfr_sfield(self, fn, None, obj, obj.fsel)
            else:
                return c.read_sfr_sfield(self, fn, self.get_cache_sbs(obj.fsel), obj)
        else:
            raise Exception("ShuffledFieldsReader.load")

    def list(self):
        return c.list_sfr(self)

    def has_sync_node(self, fn):
        return c.does_file_exist_sync_node_sfr(self, fn)

    def has(self, fn):
        if self.tags is None:
            self.tags = set(self.list())
        return fn in self.tags

class ShuffledBitSet:

    def __init__(self, fsel, new_size_node):
        assert isinstance(fsel, FieldSelection)
        self.cdata = c.mk_sbs(fsel, new_size_node)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_sbs(self)

def open_fields(path, mode, new_size_node = None):
    # path can be the folder path or the 'geon-info.txt' path
    assert isinstance(path, str)
    assert isinstance(mode, str)
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if mode == "r":
        return ShuffledFieldsReader(path, new_size_node)
    elif mode == "w":
        assert new_size_node is not None
        return ShuffledFieldsWriter(path, new_size_node)
    elif mode == "a":
        if new_size_node is None:
            return ShuffledFieldsWriter(path, [0, 0, 0, 0], True)
        else:
            return ShuffledFieldsWriter(path, new_size_node, True)
    else:
        raise Exception("open_fields")

def list_fields(path, new_size_node = None):
    assert isinstance(path, str)
    sfr = open_fields(path, "r", new_size_node)
    fns = sfr.list()
    sfr.close()
    return fns

def properly_truncate_fields(path, is_check_all = False, is_only_check = False, new_size_node = None):
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if new_size_node is None:
        return c.properly_truncate_fields_sync_node(path, is_check_all, is_only_check)
    else:
        return c.properly_truncate_fields_sync_node(path, is_check_all, is_only_check, new_size_node)

def truncate_fields(path, fns_keep, new_size_node = None):
    # fns_keep is the list of fields that need to keep
    # fns_keep needs to be in the same order as the data is stored in path
    # all fns_keep must be already in the path
    # fns_keep can be empty list
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if new_size_node is None:
        ret = c.truncate_fields_sync_node(path, fns_keep)
    else:
        ret = c.truncate_fields_sync_node(path, fns_keep, new_size_node)
    if ret != 0:
        raise Exception(f"truncate_fields error {ret}")

def check_fields(path, is_check_all = True, new_size_node = None):
    # return list of field that is stored successful
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    is_only_check = True
    return properly_truncate_fields(path, is_check_all, is_only_check, new_size_node)

def check_compressed_eigen_vectors(path):
    # return bool value suggest whether the data can be read successfully
    # return True is the data has problem
    # return False if the data is ok
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    return c.check_compressed_eigen_vectors(path)

def eigen_system_repartition(new_size_node, path, path_new = ""):
    if path[-14:] == "/geon-info.txt":
        path = path[:-14]
    if path_new[-14:] == "/geon-info.txt":
        path_new = path_new[:-14]
    return c.eigen_system_repartition(new_size_node, path, path_new)
