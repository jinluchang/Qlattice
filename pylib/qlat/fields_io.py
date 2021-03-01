import cqlat as c

from qlat.field import *
from qlat.field_selection import *

class ShuffledFieldsWriter:

    def __init__(self, path, new_size_node, is_append = False):
        assert isinstance(path, str)
        assert isinstance(is_append , bool)
        self.cdata = c.mk_sfw(path, new_size_node, is_append)
        self.fsel = None
        self.sbs = None

    def __del__(self, ):
        c.free_sfw(self)

    def new_size_node(self):
        return c.get_new_size_node_sfw(self)

    def write(self, fn, obj):
        assert isinstance(fn, str)
        if isinstance(obj, Field):
            return c.write_sfw_field(self, fn, obj)
        elif isinstance(obj, SelectedField):
            if not (self.fsel is obj.fsel):
                self.fsel = obj.fsel
                self.sbs = ShuffledBitSet(fsel, self.new_size_node())
            return c.write_sfw_sfield(self, fn, obj, self.sbs)
        else:
            raise Exception("ShuffledFieldsWriter.save")

class ShuffledFieldsReader:

    def __init__(self, path, new_size_node = None):
        assert isinstance(path, str)
        if new_size_node is None:
            self.cdata = c.mk_sfr(path)
        else:
            self.cdata = c.mk_sfr(path, new_size_node)
        # cache one fsel and its corresponding sbs
        self.fsel = None
        self.sbs = None

    def __del__(self):
        c.free_sfr(self)

    def new_size_node(self):
        return c.get_new_size_node_sfr(self)

    def read(self, fn, obj):
        assert isinstance(fn, str)
        if isinstance(obj, Field):
            return c.read_sfr_field(self, fn, obj)
        elif isinstance(obj, SelectedField):
            if not (self.fsel is obj.fsel):
                self.fsel = obj.fsel
                self.sbs = ShuffledBitSet(fsel, self.new_size_node())
            return c.read_sfr_sfield(self, fn, self.sbs, obj)
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
