cdef class Coordinate:

    def __cinit__(self):
        pass

    def __init__(self, x = None):
        if isinstance(x, Coordinate):
            self.xx = (<Coordinate>x).xx
        elif isinstance(x, (list, tuple,)):
            assert len(x) == 4
            self.xx = cc.Coordinate(x[0], x[1], x[2], x[3])
        else:
            assert x is None

    def __imatmul__(self, Coordinate v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef Coordinate x = Coordinate()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __repr__(self):
        return f"Coordinate({self.list()})"

    def list(self):
        """
        Return a list composed of the 4 components of the coordinate.
        """
        return [ self.xx[0], self.xx[1], self.xx[2], self.xx[3], ]

    def sqr(self):
        """
        Return the square sum of all the components as ``long``.
        """
        return cc.sqr(self.xx)

    def __getitem__(self, int key):
        assert 0 <= key
        assert key < 4
        return self.xx[key]

    def __setitem__(self, int key, int val):
        assert 0 <= key
        assert key < 4
        cdef int* p_val = &self.xx[key]
        p_val[0] = val
