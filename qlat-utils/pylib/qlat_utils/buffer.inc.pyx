cdef class Buffer:

    def __cinit__(self, object obj = None, int ndim = 1, int itemsize = 1):
        self.obj = obj
        self.ndim = ndim
        self.itemsize = itemsize
        self.shape_strides.resize(ndim * 2)

    cdef Py_ssize_t get_len(self):
        cdef int i
        cdef Py_ssize_t ret = 1
        for i in range(self.ndim):
            ret *= self.shape_strides[i]
        return ret

    cdef void set_strides(self):
        cdef Py_ssize_t* shapes = &self.shape_strides[0]
        cdef Py_ssize_t* strides = &self.shape_strides[self.ndim]
        cdef int i
        cdef Py_ssize_t stride = self.itemsize
        for i in range(self.ndim - 1, -1, -1):
            strides[i] = stride
            stride *= shapes[i]
