cdef class LatData:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __imatmul__(self, LatData v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef LatData x = type(self)()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef cc.bool is_complex = self.xx.is_complex()
        cdef char* fmt
        cdef int itemsize
        if is_complex:
            fmt = 'Zd'
            itemsize = sizeof(cc.Complex)
        else:
            fmt = 'd'
            itemsize = sizeof(cc.Double)
        cdef Buffer buf = Buffer(self, self.xx.ndim(), itemsize)
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(buf.ndim):
            shape[i] = self.xx.info[i].size
        buf.set_strides()
        buffer.buf = <char*>(self.xx.data())
        if flags & PyBUF_FORMAT:
            buffer.format = fmt
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = buf.itemsize
        buffer.len = buf.get_len()
        buffer.ndim = buf.ndim
        buffer.obj = buf
        buffer.readonly = 0
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

    def save_node(self, const cc.std_string& path):
        self.xx.save(path)

    def load_node(self, const cc.std_string& path):
        self.xx.load(path)

    def save(self, const cc.std_string& path):
        if cc.get_id_node() == 0:
            self.save_node(path)

    def load(self, const cc.std_string& path):
        if cc.get_id_node() == 0:
            self.load_node(path)
        self.bcast()

    def bcast(self):
        if cc.get_num_node() != 1:
            import cqlat as c
            c.bcast_lat_data(self)
        return self

    def glb_sum_in_place(self):
        if cc.get_num_node() != 1:
            import cqlat as c
            c.glb_sum_lat_data(self)
        return self

    def glb_sum(self):
        ld = self.copy()
        return ld.glb_sum_in_place()

    def __str__(self):
        return cc.show(self.xx)

    def show(self):
        return str(self)

    def set_zero(self):
        cc.set_zero(self.xx)

    def qnorm(self):
        return cc.qnorm(self.xx)

    def is_match(self, LatData ld1):
        """
        ld.info needs to be exactly equal
        """
        return cc.is_matching(self.xx, ld1.xx)

    def is_complex(self):
        return self.xx.is_complex()

    def ndim(self):
        return self.xx.ndim()

    def dim_name(self, int dim):
        assert 0 <= dim
        assert dim < self.xx.info.size()
        return self.xx.info[dim].name

    def dim_size(self, int dim):
        assert 0 <= dim
        assert dim < self.xx.info.size()
        return self.xx.info[dim].size

    def dim_indices(self, int dim):
        assert 0 <= dim
        assert dim < self.xx.info.size()
        return self.xx.info[dim].indices

    def dim_sizes(self):
        cdef int ndim = self.xx.ndim()
        cdef int i
        return [ self.xx.info[i].size for i in range(ndim) ]

    def set_dim_sizes(self, list dim_sizes, *, cc.bool is_complex = True):
        cdef int ndim = len(dim_sizes)
        cdef int ndim_real = ndim
        cdef int i
        if is_complex:
            ndim_real += 1
        self.xx.info.resize(ndim_real)
        for i in range(ndim):
            self.xx.info[i].size = dim_sizes[i]
        if is_complex:
            self.xx.info[ndim] = cc.lat_dim_re_im()
        cc.lat_data_alloc(self.xx)

    def set_dim_name(self, int dim, const cc.std_string& name, list indices = None):
        assert 0 <= dim
        assert dim < self.xx.info.size()
        cdef int size
        cdef int i
        self.xx.info[dim].name = name
        if indices is None:
            self.xx.info[dim].indices.resize(0)
        else:
            indices = [ str(idx).replace("\n", "  ") for idx in indices ]
            size = len(indices)
            self.xx.info[dim].indices.resize(size)
            for i in range(size):
                self.xx.info[dim].indices[i] = indices[i]

    def dim_names(self):
        """
        by default, return list can be used as the input argument for ld.from_numpy
        """
        cdef int ndim = self.xx.ndim()
        cdef int i
        return [ self.xx.info[i].name for i in range(ndim) ]

    def to_numpy(self):
        return np.asarray(self).copy()

    def from_numpy(self, numpy.ndarray val, list dim_names = None, *, cc.bool is_complex = True):
        """
        only set LatData shape if it is initially empty
        otherwise only set data and ignore shape completely
        dim_names should be a list of names for each dimension
        """
        cdef int ndim = val.ndim
        cdef int dim
        cdef list shape = [ val.shape[dim] for dim in range(ndim) ]
        if self.ndim() == 0:
            self.set_dim_sizes(shape, is_complex = is_complex)
            if dim_names is None:
                dim_names = [ n for n in "ijklmnopqrstuvwxyz" ]
            assert ndim <= len(dim_names)
            for dim in range(ndim):
                name = dim_names[dim]
                self.set_dim_name(dim, name)
        np.asarray(self).ravel()[:] = val.ravel()[:]
        return self

    def to_list(self):
        return np.asarray(self).ravel().tolist()

    def from_list(self, list val, *, cc.bool is_complex = True):
        if self.ndim() == 0:
            self.set_dim_sizes([ len(val), ], is_complex = is_complex)
            self.set_dim_name(0, "i")
        np.asarray(self).ravel()[:] = np.array(val)
        return self

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

    def info(self, dim = None):
        """
        by default, return list can be used as the input argument for ld.set_info or mk_lat_data
        """
        if dim is None:
            ndim = self.ndim()
            return [ self.info(i) for i in range(ndim) ]
        else:
            dim_name = self.dim_name(dim)
            dim_size = self.dim_size(dim)
            dim_indices = self.dim_indices(dim)
            return [ dim_name, dim_size, dim_indices, ]

    def set_info(self, list info_list, *, cc.bool is_complex = True):
        """
        ``info_list`` format::\n
            [ [ dim_name, dim_size, dim_indices, ], ... ]
        dim_indices can be optional
        """
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
                raise Exception(f"LatData setinfo info_list={info_list}")
        self.set_zero()

    def __iadd__(self, LatData ld1):
        self.xx = self.xx + ld1.xx
        return self

    def __isub__(self, LatData ld1):
        self.xx = self.xx - ld1.xx
        return self

    def __imul__(self, factor):
        cdef cc.Double f
        cdef cc.Complex c
        if isinstance(factor, float):
            f = factor
            self.xx = self.xx * f
        elif isinstance(factor, complex):
            c = factor
            self.xx = self.xx * c
        else:
            assert False
        return self

    def __add__(LatData ld1, LatData ld2):
        cdef LatData ld = type(ld1)()
        ld.xx = ld1.xx + ld2.xx
        return ld

    def __sub__(LatData ld1, LatData ld2):
        cdef LatData ld = type(ld1)()
        ld.xx = ld1.xx - ld2.xx
        return ld

    def __mul__(ld, factor):
        if isinstance(factor, LatData):
            ld, factor = factor, ld
        assert isinstance(ld, LatData)
        cdef LatData ld1 = type(ld)()
        cdef LatData ld0 = ld
        cdef cc.Double f
        cdef cc.Complex c
        if isinstance(factor, float):
            f = factor
            ld1.xx = ld0.xx * f
        elif isinstance(factor, complex):
            c = factor
            ld1.xx = ld0.xx * c
        else:
            assert False
        return ld1

    def __neg__(self):
        cdef LatData ld = type(self)()
        ld.xx = ld.xx - self.xx
        return ld

    def __pos__(self):
        return self

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

### -------------------------------------------------------------------

def mk_lat_data(list info_list, *, cc.bool is_complex = True):
    """
    ``info_list`` format::\n
        [ [ dim_name, dim_size, dim_indices, ], ... ]
    dim_indices can be optional
    """
    ld = LatData()
    ld.set_info(info_list, is_complex = is_complex)
    return ld

def load_lat_data(const cc.std_string& path):
    """
    Load ``lat_data`` from file ``path``.
    """
    ld = LatData()
    ld.load(path)
    return ld
