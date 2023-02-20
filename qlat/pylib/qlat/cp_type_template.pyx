### -------------------------------------------------------------------

cdef class FieldTYPENAME(FieldBase):

    ctype = ElemTypeTYPENAME

    def __cinit__(self):
        self.cdata = <long>&(self.xx)
        self.view_count = 0

    def __init__(self, Geometry geo = None, int multiplicity = 0):
        if geo is None:
            return
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        self.xx.init(geo.xx, multiplicity)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef const cc.Geometry* p_geo = &self.xx.get_geo()
        cdef cc.Coordinate local_site = p_geo[0].local_site()
        cdef int multiplicity = p_geo[0].multiplicity
        cdef cc.Vector[cc.TYPENAME] fvec = cc.get_data(self.xx)
        cdef int ndim = 5 + ElemTypeTYPENAME.ndim()
        cdef char* fmt = ElemTypeTYPENAME.format()
        cdef Buffer buf = Buffer(self, ndim, ElemTypeTYPENAME.itemsize())
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeTYPENAME.shape()
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(4):
            shape[i] = local_site[i]
        shape[4] = multiplicity
        for i in range(ElemTypeTYPENAME.ndim()):
            shape[5 + i] = vec[i]
        buf.set_strides()
        buffer.buf = <char*>(fvec.data())
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
        assert buffer.len * buffer.itemsize == fvec.size() * ElemTypeTYPENAME.size()
        self.view_count += 1

    def __releasebuffer__(self, Py_buffer *buffer):
        self.view_count -= 1

    def __imatmul__(self, f1):
        # f1 can be Field, SelectedField, SelectedPoints
        # field geo does not change if already initialized
        if isinstance(f1, FieldTYPENAME):
            self.xx = (<FieldTYPENAME>f1).xx
        elif isinstance(f1, SelectedFieldTYPENAME):
            c.set_field_sfield(self, f1)
        elif isinstance(f1, SelectedPointsTYPENAME):
            c.set_field_spfield(self, f1)
        else:
            raise Exception(f"FieldTYPENAME @= type mismatch {type(self)} {type(f1)}")
        return self

    @q.timer
    def copy(self, is_copying_data = True):
        f = type(self)()
        if is_copying_data:
            (<FieldTYPENAME>f).xx = self.xx
        return f

    @q.timer
    def set_zero(self):
        cc.set_zero(self.xx)

    def swap(self, FieldTYPENAME f1):
        cc.qswap(f1.xx, self.xx)

### -------------------------------------------------------------------

cdef class SelectedFieldTYPENAME(SelectedFieldBase):

    ctype = ElemTypeTYPENAME

    def __cinit__(self):
        self.cdata = <long>&(self.xx)
        self.view_count = 0

    def __init__(self, FieldSelection fsel, int multiplicity = 0):
        self.fsel = fsel
        if multiplicity > 0 and self.fsel is not None:
            if self.view_count > 0:
                raise ValueError("can't re-init while being viewed")
            self.xx.init(self.fsel.xx, multiplicity)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef long n_elems = self.xx.n_elems
        cdef int multiplicity = self.xx.get_geo().multiplicity
        cdef cc.Vector[cc.TYPENAME] fvec = cc.get_data(self.xx)
        cdef int ndim = 2 + ElemTypeTYPENAME.ndim()
        cdef char* fmt = ElemTypeTYPENAME.format()
        cdef Buffer buf = Buffer(self, ndim, ElemTypeTYPENAME.itemsize())
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeTYPENAME.shape()
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        shape[0] = n_elems
        shape[1] = multiplicity
        for i in range(ElemTypeTYPENAME.ndim()):
            shape[2 + i] = vec[i]
        buf.set_strides()
        buffer.buf = <char*>(fvec.data())
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
        assert buffer.len * buffer.itemsize == fvec.size() * ElemTypeTYPENAME.size()
        self.view_count += 1

    def __releasebuffer__(self, Py_buffer *buffer):
        self.view_count -= 1

    def __imatmul__(self, f1):
        # won't change self.fsel
        if isinstance(f1, SelectedFieldTYPENAME):
            # two fsel do not need to match
            if self.fsel is f1.fsel:
                self.xx = (<SelectedFieldTYPENAME>f1).xx
            else:
                c.set_sfield_sfield(self, f1)
        elif isinstance(f1, FieldTYPENAME):
            c.set_sfield_field(self, f1)
        elif isinstance(f1, SelectedPointsTYPENAME):
            # psel may not have to be subset of fsel
            c.set_sfield_spfield(self, f1)
        else:
            raise Exception(f"SelectedFieldTYPENAME @= type mismatch {type(self)} {type(f1)}")
        return self

    @q.timer
    def copy(self, is_copying_data = True):
        f = type(self)(self.fsel)
        if is_copying_data:
            (<SelectedFieldTYPENAME>f).xx = self.xx
        return f

    @q.timer
    def set_zero(self):
        cc.set_zero(self.xx)

    def swap(self, SelectedFieldTYPENAME f1):
        assert f1.fsel is self.fsel
        cc.qswap(f1.xx, self.xx)

### -------------------------------------------------------------------

cdef class SelectedPointsTYPENAME(SelectedPointsBase):

    ctype = ElemTypeTYPENAME

    def __cinit__(self):
        self.cdata = <long>&(self.xx)
        self.view_count = 0

    def __init__(self, PointSelection psel, int multiplicity = 0):
        self.psel = psel
        if multiplicity > 0 and self.psel is not None:
            if self.view_count > 0:
                raise ValueError("can't re-init while being viewed")
            self.xx.init(self.psel.xx, multiplicity)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef long n_points = self.xx.n_points
        cdef int multiplicity = self.xx.multiplicity
        cdef cc.Vector[cc.TYPENAME] fvec = cc.get_data(self.xx)
        cdef int ndim = 2 + ElemTypeTYPENAME.ndim()
        cdef char* fmt = ElemTypeTYPENAME.format()
        cdef Buffer buf = Buffer(self, ndim, ElemTypeTYPENAME.itemsize())
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeTYPENAME.shape()
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        shape[0] = n_points
        shape[1] = multiplicity
        for i in range(ElemTypeTYPENAME.ndim()):
            shape[2 + i] = vec[i]
        buf.set_strides()
        buffer.buf = <char*>(fvec.data())
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
        assert buffer.len * buffer.itemsize == fvec.size() * ElemTypeTYPENAME.size()
        self.view_count += 1

    def __releasebuffer__(self, Py_buffer *buffer):
        self.view_count -= 1

    def __imatmul__(self, f1):
        # won't change self.psel
        if isinstance(f1, SelectedPointsTYPENAME):
            # two psel must be the same object
            if self.psel is f1.psel:
                self.xx = (<SelectedPointsTYPENAME>f1).xx
            else:
                raise Exception("SelectedPointsTYPENAME @= psel mismatch")
        elif isinstance(f1, SelectedFieldTYPENAME):
            # only assign available points
            c.set_spfield_sfield(self, f1)
        elif isinstance(f1, FieldTYPENAME):
            c.set_spfield_field(self, f1)
        else:
            raise Exception(f"SelectedPointsTYPENAME @= type mismatch {type(self)} {type(f1)}")
        return self

    @q.timer
    def copy(self, is_copying_data = True):
        f = type(self)(self.psel)
        if is_copying_data:
            (<SelectedPointsTYPENAME>f).xx = self.xx
        return f

    @q.timer
    def set_zero(self):
        cc.set_zero(self.xx)

    def swap(self, SelectedPointsTYPENAME f1):
        assert f1.psel is self.psel
        cc.qswap(f1.xx, self.xx)

### -------------------------------------------------------------------

field_type_dict[ElemTypeTYPENAME] = FieldTYPENAME

selected_field_type_dict[ElemTypeTYPENAME] = SelectedFieldTYPENAME

selected_points_type_dict[ElemTypeTYPENAME] = SelectedPointsTYPENAME

### -------------------------------------------------------------------
