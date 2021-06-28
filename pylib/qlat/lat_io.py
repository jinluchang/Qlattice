import cqlat as c

class LatData:

    def __init__(self):
        self.cdata = c.mk_lat_data()

    def __del__(self):
        c.free_lat_data(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, LatData)
        c.set_lat_data(self, v1)
        return self

    def copy(self):
        x = LatData()
        x @= self
        return x

    def save_node(self, path):
        c.save_lat_data(self, path)

    def load_node(self, path):
        c.load_lat_data(self, path)

    def bcast(self):
        c.bcast_lat_data(self)

    def save(self, path):
        from qlat.mpi import get_id_node
        if get_id_node() == 0:
            self.save_node(path)

    def load(self, path):
        from qlat.mpi import get_id_node
        if get_id_node() == 0:
            self.load_node(path)
        self.bcast()

    def show(self):
        return c.show_lat_data(self)

    def __iadd__(self, ld1):
        assert isinstance(ld1, LatData)
        c.set_add_lat_data(self, ld1)
        return self

    def __isub__(self, ld1):
        assert isinstance(ld1, LatData)
        c.set_sub_lat_data(self, ld1)
        return self

    def __imul__(self, factor):
        assert isinstance(factor, float)
        c.set_mul_double_lat_data(self, factor)
        return self

    def set_zero(self):
        return c.set_zero_lat_data(self)

    def qnorm(self):
        return c.qnorm_lat_data(self)

    def is_match(self, ld1):
        # ld.info needs to be exactly equal
        return c.is_matching_lat_data(self, ld1)

    def is_complex(self):
        return c.is_complex_lat_data(self)

    def ndim(self):
        return c.get_ndim_lat_data(self)

    def dim_sizes(self):
        return c.get_dim_sizes_lat_data(self)

    def dim_name(self, dim):
        return c.get_dim_name_lat_data(self, dim)

    def dim_indices(self, dim):
        return c.get_dim_indices_lat_data(self, dim)

    def set_dim_sizes(self, dim_sizes, *, is_complex = True):
        return c.set_dim_sizes_lat_data(self, dim_sizes, is_complex)

    def set_dim_name(self, dim, name, indices = []):
        return c.set_dim_name_lat_data(self, dim, name, indices)

    def to_list(self, *, is_complex = True):
        is_always_double = not is_complex
        return c.peek_lat_data(self, [], is_always_double)

    def from_list(self, val, *, is_complex = True):
        if self.ndim() == 0:
            self.set_dim_sizes([len(val)], is_complex)
            self.set_dim_name(0, "i")
        is_always_double = not is_complex
        c.poke_lat_data(self, [], val, is_always_double)
        return self

    def __setitem__(self, idx, val):
        # use list with correct length as val
        return c.poke_lat_data(self, idx, val)

    def __getitem__(self, idx):
        # return a new list every call
        return c.peek_lat_data(self, idx)

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
