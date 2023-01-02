### -------------------------------------------------------------------

cdef class FieldTYPENAME(FieldBase):

    ctype = ElemTypeTYPENAME

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __init__(self, geo = None, multiplicity = None):
        if geo is None:
            return
        assert isinstance(geo, Geometry)
        if multiplicity is None:
            self.xx.init((<Geometry>geo).xx)
        else:
            assert isinstance(multiplicity, int)
            self.xx.init((<Geometry>geo).xx, <int>multiplicity)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef cc.Coordinate local_site = self.xx.get_geo().local_site()
        cdef int ndim = 4 + ElemTypeTYPENAME.ndim()
        cdef char* fmt = ElemTypeTYPENAME.format()
        cdef Buffer buf = Buffer(self, ndim, ElemTypeTYPENAME.itemsize())
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeTYPENAME.shape()
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(4):
            shape[i] = local_site[i]
        for i in range(buf.ndim):
            shape[4 + i] = vec[i]
        buf.set_strides()
        buffer.buf = <char*>(cc.get_data(self.xx).data())
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

    def __imatmul__(self, f1):
        # f1 can be Field, SelectedField, SelectedPoints
        # field geo does not change if already initialized
        if isinstance(f1, FieldTYPENAME):
            self.xx = f1.xx
        else:
            assert f1.ctype is self.ctype
            from qlat.selected_field import SelectedField
            from qlat.selected_points import SelectedPoints
            if isinstance(f1, SelectedField):
                c.set_field_sfield(self, f1)
            elif isinstance(f1, SelectedPoints):
                c.set_field_spfield(self, f1)
            else:
                raise Exception(f"Field @= type mismatch {type(self)} {type(f1)}")
        return self

    def copy(self, is_copying_data = True):
        f = FieldTYPENAME()
        if is_copying_data:
            f @= self
        return f

    def swap(self, FieldTYPENAME f1):
        cc.qswap(f1.xx, self.xx)

### -------------------------------------------------------------------

field_type_dict[ElemTypeTYPENAME] = FieldTYPENAME

### -------------------------------------------------------------------
