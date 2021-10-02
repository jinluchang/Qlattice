import cqlat as c

from qlat.field import *
from qlat.cache import *
from qlat.field_selection import *
from qlat.selected_field import *
from qlat.utils_io import *

cache_fields_io = mk_cache("fields_io")

class ShuffledFieldsWriter:

    def __init__(self, path, new_size_node, is_append = False):
        assert isinstance(path, str)
        assert isinstance(is_append , bool)
        mk_file_dirs_info(path)
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
        if isinstance(obj, Field):
            return c.write_sfw_field(self, fn, obj)
        elif isinstance(obj, SelectedField):
            return c.write_sfw_sfield(self, fn, obj, self.get_cache_sbs(obj.fsel))
        else:
            raise Exception("ShuffledFieldsWriter.save")

    def flush(self):
        return c.flush_sfw(self)

class ShuffledFieldsReader:

    def __init__(self, path, new_size_node = None):
        assert isinstance(path, str)
        if new_size_node is None:
            self.cdata = c.mk_sfr(path)
        else:
            self.cdata = c.mk_sfr(path, new_size_node)

    def close(self):
        if self.cdata is not None:
            c.free_sfr(self)
        self.cdata = None
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
        assert isinstance(fn, str)
        if isinstance(obj, Field):
            return c.read_sfr_field(self, fn, obj)
        elif isinstance(obj, SelectedField):
            return c.read_sfr_sfield(self, fn, self.get_cache_sbs(obj.fsel), obj)
        else:
            raise Exception("ShuffledFieldsReader.load")

    def list(self):
        return c.list_sfr(self)

class ShuffledBitSet:

    def __init__(self, fsel, new_size_node):
        assert isinstance(fsel, FieldSelection)
        self.cdata = c.mk_sbs(fsel, new_size_node)

    def __del__(self):
        c.free_sbs(self)

def open_fields(path, mode, new_size_node = None):
    assert isinstance(path, str)
    assert isinstance(mode, str)
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
    sfr = open_fields(path, "r", new_size_node)
    fns = sfr.list()
    sfr.close()
    return fns

def properly_truncate_fields(path, is_check_all = False, is_only_check = False, new_size_node = None):
    if new_size_node is None:
        return c.properly_truncate_fields_sync_node(path, is_check_all, is_only_check)
    else:
        return c.properly_truncate_fields_sync_node(path, is_check_all, is_only_check, new_size_node)

def check_fields(path, is_check_all = True, new_size_node = None):
    # return list of field that is stored successful
    is_only_check = True
    return properly_truncate_fields(path, is_check_all, is_only_check, new_size_node)

def check_compressed_eigen_vectors(path):
    # return bool value suggest whether the data can be read successfully
    return c.check_compressed_eigen_vectors(path)
