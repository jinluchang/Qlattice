import qlat.cqlat as c

from qlat_utils import *
from qlat.geometry import *
from qlat.field import *
from qlat.field_selection import *
from qlat.utils_io import *

class SelectedField:

    # self.ctype
    # self.fsel
    # self.cdata

    def __init__(self, ctype, fsel, multiplicity = None):
        assert isinstance(ctype, str)
        assert fsel is None or isinstance(fsel, FieldSelection)
        self.ctype = ctype
        self.fsel = fsel
        if fsel is None or multiplicity is None:
            self.cdata = c.mk_sfield(ctype)
        else:
            assert isinstance(multiplicity, int)
            self.cdata = c.mk_sfield_fsel(ctype, fsel, multiplicity)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_sfield(self)

    def __imatmul__(self, f1):
        # won't change self.fsel
        assert f1.ctype == self.ctype
        if isinstance(f1, SelectedField):
            # two fsel do not need to match
            if self.fsel is f1.fsel:
                c.set_sfield(self, f1)
            else:
                c.set_sfield_sfield(self, f1)
        elif isinstance(f1, Field):
            c.set_sfield_field(self, f1)
        elif isinstance(f1, SelectedPoints):
            # psel must be subset of fsel
            c.set_sfield_spfield(self, f1)
        else:
            assert False
        return self

    def copy(self, is_copying_data = True):
        f = SelectedField(self.ctype, self.fsel)
        if is_copying_data:
            f @= self
        return f

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def swap(self, x):
        assert isinstance(x, SelectedField)
        assert x.ctype == self.ctype
        assert x.fsel is self.fsel
        cdata = x.cdata
        x.cdata = self.cdata
        self.cdata = cdata

    def n_elems(self):
        return c.get_n_elems_sfield(self)

    def total_site(self):
        return c.get_total_site_sfield(self)

    def multiplicity(self):
        return c.get_multiplicity_sfield(self)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
        c.set_geo_sfield(geo, self)
        return geo

    def __iadd__(self, f1):
        assert isinstance(f1, SelectedField)
        assert f1.ctype == self.ctype
        c.set_add_sfield(self, f1)
        return self

    def __isub__(self, f1):
        assert isinstance(f1, SelectedField)
        assert f1.ctype == self.ctype
        c.set_sub_sfield(self, f1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_sfield(self, factor)
        return self

    def set_zero(self):
        c.set_zero_sfield(self)

    def qnorm(self):
        return c.qnorm_sfield(self)

    def get_elems(self, idx):
        return np.array(c.get_elems_sfield(self, idx))

    def get_elem(self, idx, m = None):
        if m is None:
            return np.array(c.get_elem_sfield(self, idx))
        else:
            return np.array(c.get_elem_sfield(self, idx, m))

    def set_elems(self, idx, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elems_sfield(self, idx, val)
        elif isinstance(val, np.ndarray):
            return self.set_elems(idx, val.tobytes())
        else:
            assert False

    def set_elem(self, idx, m, val):
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(val, bytes):
            return c.set_elem_sfield(self, idx, m, val)
        elif isinstance(val, np.ndarray):
            return self.set_elem(idx, m, val.tobytes())
        else:
            assert False

    def __getitem__(self, i):
        # i can be (idx, m,) or idx
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], int):
            idx, m = i
            return self.get_elem(idx, m)
        elif isinstance(i, int):
            idx = i
            return self.get_elems(idx)
        else:
            assert False
            return None

    def __setitem__(self, i, val):
        # i can be (idx, m,) or idx
        # val should be np.ndarray or bytes. e.g. np.array([1, 2, 3], dtype = complex).tobytes()
        if isinstance(i, tuple) and len(i) == 2 and isinstance(i[0], int):
            idx, m = i
            return self.set_elem(idx, m, val)
        elif isinstance(i, int):
            idx = i
            return self.set_elems(idx, val)
        else:
            assert False
            return None

    def sparse(self, sel):
        # deprecated
        if isinstance(sel, PointSelection):
            from qlat.selected_points import SelectedPoints, set_selected_points
            psel = sel
            sp = SelectedPoints(self.ctype, psel)
            set_selected_points(sp, self)
            return sp
        elif isinstance(sel, FieldSelection):
            fsel = sel
            sf = SelectedField(self.ctype, fsel)
            set_selected_field(sf, self)
            return sf
        else:
            raise Exception("Field.sparse")

    def save(self, path, *args):
        # possible way to call:
        # f.save(path)
        # f.save(path, new_size_node)
        # f.save(sfw, fn)
        from qlat.fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            assert len(args) == 0
            mk_file_dirs_info(path)
            return c.save_sfield(self, path)
        elif isinstance(path, ShuffledFieldsWriter):
            sfw = path
            [fn] = args
            return sfw.write(fn, self)
        else:
            raise Exception("SelectedField.save")

    def load(self, path, *args):
        # possible way to call:
        # f.load(path)
        # f.load(sfr, fn)
        # if self.fsel is None, self.fsel will be set during f.load(sfr, fn)
        from qlat.fields_io import ShuffledFieldsReader
        if isinstance(path, str):
            return c.load_sfield(self, path)
        elif isinstance(path, ShuffledFieldsReader):
            sfr = path
            [fn] = args
            return sfr.read(fn, self)
        else:
            raise Exception("SelectedField.load")

    def save_64(self, path, *args):
        f = self.copy()
        from qlat.fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            f.to_from_endianness("big_64")
        elif isinstance(path, ShuffledFieldsWriter):
            f.to_from_endianness("little_64")
        return f.save(path, *args)

    def save_double(self, path, *args):
        return self.save_64(path, *args)

    def save_float_from_double(self, path, *args):
        ff = SelectedField("float", self.fsel)
        ff.float_from_double(self)
        from qlat.fields_io import ShuffledFieldsWriter
        if isinstance(path, str):
            ff.to_from_endianness("big_32")
        elif isinstance(path, ShuffledFieldsWriter):
            ff.to_from_endianness("little_32")
        return ff.save(path, *args)

    def load_64(self, path, *args):
        ret = self.load(path, *args)
        if ret > 0:
            from qlat.fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                self.to_from_endianness("big_64")
            elif isinstance(path, ShuffledFieldsReader):
                self.to_from_endianness("little_64")
        return ret

    def load_double(self, path, *args):
        return self.load_64(path, *args)

    def load_double_from_float(self, path, *args):
        ff = SelectedField("float", self.fsel)
        ret = ff.load(path, *args)
        if ret > 0:
            from qlat.fields_io import ShuffledFieldsReader
            if isinstance(path, str):
                ff.to_from_endianness("big_32")
            elif isinstance(path, ShuffledFieldsReader):
                ff.to_from_endianness("little_32")
            self.double_from_float(ff)
        return ret

    def float_from_double(self, f):
        assert isinstance(f, SelectedField)
        assert self.ctype == "float"
        self.fsel = f.fsel
        c.convert_float_from_double_sfield(self, f)

    def double_from_float(self, ff):
        assert isinstance(ff, SelectedField)
        assert ff.ctype == "float"
        self.fsel = ff.fsel
        c.convert_double_from_float_sfield(self, ff)

    def to_from_endianness(self, tag):
        assert isinstance(tag, str)
        c.to_from_endianness_sfield(self, tag)

    def field_shift(self, shift, is_reflect = False):
        # return new shifted SelectedField
        # shift is the coordinate to shift the field
        # is_reflect determine whether to negate coordinate after shift
        f1 = self.copy(is_copying_data = False)
        f1.fsel = FieldSelection()
        c.field_shift_sfield(f1, self, shift, is_reflect)
        return f1

    def glb_sum_tslice(self, *, t_dir = 3):
        # return SelectedPoints(self.ctype, get_psel_tslice(self.total_site(), t_dir = t_dir))
        from qlat.selected_points import SelectedPoints
        psel = get_psel_tslice(self.total_site(), t_dir = t_dir)
        sp = SelectedPoints(self.ctype, psel)
        if self.ctype in field_ctypes_double:
            c.glb_sum_tslice_double_sfield(sp, self, t_dir)
        elif self.ctype in field_ctypes_long:
            c.glb_sum_tslice_long_sfield(sp, self, t_dir)
        else:
            assert False
        return sp

###

@timer
def set_selected_field(sf, f):
    # deprecated use @=
    assert isinstance(sf, SelectedField)
    if isinstance(f, Field):
        c.set_sfield_field(sf, f)
    elif isinstance(f, SelectedField):
        c.set_sfield_sfield(sf, f)
    else:
        raise Exception("set_selected_field")
