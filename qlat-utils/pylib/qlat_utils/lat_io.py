import cqlat_utils as cu

import numpy as np

from qlat_utils.utils_io import *

class LatData:

    # self.cdata

    def __init__(self):
        self.cdata = cu.mk_lat_data()

    def __del__(self):
        assert isinstance(self.cdata, int)
        cu.free_lat_data(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, LatData)
        cu.set_lat_data(self, v1)
        return self

    def copy(self, is_copying_data = True):
        x = LatData()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def save_node(self, path):
        mk_file_dirs(path)
        cu.save_lat_data(self, path)

    def load_node(self, path):
        cu.load_lat_data(self, path)

    def bcast(self):
        if get_num_node() != 1:
            import cqlat as c
            c.bcast_lat_data(self)
        return self

    def glb_sum_in_place(self):
        if get_num_node() != 1:
            import cqlat as c
            c.glb_sum_lat_data(self)
        return self

    def glb_sum(self):
        ld = self.copy()
        return ld.glb_sum_in_place()

    def save(self, path):
        if get_id_node() == 0:
            mk_file_dirs(path)
            self.save_node(path)

    def load(self, path):
        if get_id_node() == 0:
            self.load_node(path)
        self.bcast()

    def show(self):
        return cu.show_lat_data(self)

    def __iadd__(self, ld1):
        assert isinstance(ld1, LatData)
        cu.set_add_lat_data(self, ld1)
        return self

    def __isub__(self, ld1):
        assert isinstance(ld1, LatData)
        cu.set_sub_lat_data(self, ld1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        cu.set_mul_double_lat_data(self, factor)
        return self

    def __add__(self, ld1):
        ld = self.copy()
        ld += ld1
        return ld

    def __sub__(self, ld1):
        ld = self.copy()
        ld -= ld1
        return ld

    def __mul__(self, factor):
        ld = self.copy()
        ld *= factor
        return ld

    __rmul__ = __mul__

    def __neg__(self):
        ld = self.copy()
        ld *= -1
        return ld

    def __pos__(self):
        return self

    def set_zero(self):
        return cu.set_zero_lat_data(self)

    def qnorm(self):
        return cu.qnorm_lat_data(self)

    def is_match(self, ld1):
        # ld.info needs to be exactly equal
        return cu.is_matching_lat_data(self, ld1)

    def is_complex(self):
        return cu.is_complex_lat_data(self)

    def ndim(self, *, is_complex = True):
        is_always_double = not is_complex
        return cu.get_ndim_lat_data(self, is_always_double)

    def dim_sizes(self, *, is_complex = True):
        is_always_double = not is_complex
        return cu.get_dim_sizes_lat_data(self, is_always_double)

    def dim_name(self, dim):
        return cu.get_dim_name_lat_data(self, dim)

    def dim_size(self, dim):
        return cu.get_dim_size_lat_data(self, dim)

    def dim_indices(self, dim):
        return cu.get_dim_indices_lat_data(self, dim)

    def set_dim_sizes(self, dim_sizes, *, is_complex = True):
        return cu.set_dim_sizes_lat_data(self, dim_sizes, is_complex)

    def set_dim_name(self, dim, name, indices = None):
        if indices is None:
            indices = []
        else:
            indices = [ str(idx).replace("\n", "  ") for idx in indices ]
        return cu.set_dim_name_lat_data(self, dim, name, indices)

    def dim_names(self, *, is_complex = True):
        # by default, return list can be used as the input argument for ld.from_numpy
        ndim = self.ndim(is_complex = is_complex)
        return [ self.dim_name(dim) for dim in range(ndim) ]

    def to_list(self, *, is_complex = True):
        is_always_double = not is_complex
        return cu.peek_lat_data(self, [], is_always_double)

    def from_list(self, val, *, is_complex = True):
        if self.ndim() == 0:
            self.set_dim_sizes([len(val)], is_complex = is_complex)
            self.set_dim_name(0, "i")
        is_always_double = not is_complex
        cu.poke_lat_data(self, [], val, is_always_double)
        return self

    def to_numpy(self, *, is_complex = True):
        is_always_double = not is_complex
        v = np.array(cu.peek_lat_data(self, [], is_always_double))
        return v.reshape(self.dim_sizes(is_complex = is_complex))

    def from_numpy(self, val, dim_names = None, *, is_complex = True):
        # only set LatData shape if it is initially empty
        # otherwise only set data and ignore shape completely
        # dim_names should be a list of names for each dimension
        if self.ndim() == 0:
            if dim_names is None:
                dim_names = "ijklmnopqrstuvwxyz"
            self.set_dim_sizes(list(val.shape), is_complex = is_complex)
            for dim, (dummy_size, name) in enumerate(zip(val.shape, dim_names)):
                self.set_dim_name(dim, name)
        is_always_double = not is_complex
        cu.poke_lat_data(self, [], list(val.flatten()), is_always_double)
        return self

    def __setitem__(self, idx, val):
        # use list with correct length as val
        # idx should be tuple or list of int
        if isinstance(idx, int):
            idx = [ idx, ]
        if isinstance(val, np.ndarray):
            val = val.flatten()
        return cu.poke_lat_data(self, idx, list(val))

    def __getitem__(self, idx):
        # return a new list every call
        # idx should be tuple or list of int
        if isinstance(idx, int):
            idx = [ idx, ]
        shape = self.dim_sizes()[len(idx):]
        return np.array(cu.peek_lat_data(self, idx)).reshape(shape)

    def __getstate__(self):
        is_complex = self.is_complex()
        ndim = self.ndim()
        dim_sizes = self.dim_sizes()
        assert len(dim_sizes) == ndim
        dim_names = [ self.dim_name(dim) for dim in range(ndim) ]
        dim_indices = [ self.dim_indices(dim) for dim in range(ndim) ]
        data_list = self.to_list()
        return [ is_complex, dim_sizes, dim_names, dim_indices, data_list ]

    def __setstate__(self, state):
        [ is_complex, dim_sizes, dim_names, dim_indices, data_list ] = state
        self.__init__()
        self.set_dim_sizes(dim_sizes, is_complex = is_complex)
        ndim = len(dim_sizes)
        for dim in range(ndim):
            self.set_dim_name(dim, dim_names[dim], dim_indices[dim])
        self.from_list(data_list)

    def info(self, dim = None, *, is_complex = True):
        # by default, return list can be used as the input argument for ld.set_info or mk_lat_data
        if dim is None:
            ndim = self.ndim(is_complex = is_complex)
            return [ self.info(i) for i in range(ndim) ]
        else:
            dim_name = self.dim_name(dim)
            dim_size = self.dim_size(dim)
            dim_indices = self.dim_indices(dim)
            return [ dim_name, dim_size, dim_indices, ]

    def set_info(self, info_list, *, is_complex = True):
        # info_list format:
        # [ [ dim_name, dim_size, dim_indices, ], ... ]
        # dim_indices can be optional
        for info in info_list:
            assert len(info) >= 2
        dim_sizes = [ info[1] for info in info_list ]
        self.set_dim_sizes(dim_sizes, is_complex = is_complex)
        ndim = len(dim_sizes)
        for dim in range(ndim):
            info = info_list[dim]
            if len(info) == 2:
                dim_name, dummy_dim_size = info
                self.set_dim_name(dim, dim_name)
            elif len(info) == 3:
                dim_name, dummy_dim_size, dim_indices = info
                self.set_dim_name(dim, dim_name, dim_indices)
            else:
                raise Exception(f"LatData setinfo info_list={info_list} is_complex={is_complex}")
        self.set_zero()

###

def mk_lat_data(info_list, *, is_complex = True):
    ld = LatData()
    ld.set_info(info_list, is_complex = is_complex)
    return ld

def load_lat_data(path):
    ld = LatData()
    ld.load(path)
    return ld
