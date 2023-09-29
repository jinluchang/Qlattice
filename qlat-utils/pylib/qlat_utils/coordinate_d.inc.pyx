cdef class CoordinateD:

    def __cinit__(self):
        pass

    def __init__(self, x = None):
        if isinstance(x, CoordinateD):
            self.xx = (<CoordinateD>x).xx
        elif isinstance(x, Coordinate):
            self.xx = cc.CoordinateD((<Coordinate>x).xx)
        elif isinstance(x, (list, tuple,)):
            assert len(x) == 4
            self.xx = cc.CoordinateD(x[0], x[1], x[2], x[3])
        else:
            assert x is None

    def __imatmul__(self, CoordinateD v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef CoordinateD x = CoordinateD()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __repr__(self):
        return f"CoordinateD({self.to_list()})"

    def to_list(self):
        """
        Return a list composed of the 4 components of the coordinate.
        """
        return [ self.xx[0], self.xx[1], self.xx[2], self.xx[3], ]

    def from_list(self, list x):
        """
        set value based on a list composed of the 4 components of the coordinate.
        """
        assert isinstance(x, (list, tuple,))
        assert len(x) == 4
        self.xx = cc.CoordinateD(x[0], x[1], x[2], x[3])

    def __getstate__(self):
        return self.to_list()

    def __setstate__(self, state):
        self.from_list(state)

    def __getitem__(self, int key):
        assert 0 <= key
        assert key < 4
        return self.xx[key]

    def __setitem__(self, int key, double val):
        assert 0 <= key
        assert key < 4
        cdef double* p_val = &self.xx[key]
        p_val[0] = val

    def __add__(CoordinateD c1, CoordinateD c2):
        cdef CoordinateD x = CoordinateD()
        x.xx = c1.xx + c2.xx
        return x

    def __sub__(CoordinateD c1, CoordinateD c2):
        cdef CoordinateD x = CoordinateD()
        x.xx = c1.xx - c2.xx
        return x

    def __mul__(c1, c2):
        cdef CoordinateD x = CoordinateD()
        if isinstance(c1, CoordinateD) and isinstance(c2, (int, float,)):
            x.xx = (<CoordinateD>c1).xx * (<double>c2)
        elif isinstance(c1, (int, float,)) and isinstance(c2, CoordinateD):
            x.xx = (<double>c1) * (<CoordinateD>c2).xx
        elif isinstance(c1, CoordinateD) and isinstance(c2, CoordinateD):
            x.xx = (<CoordinateD>c1).xx * (<CoordinateD>c2).xx
        else:
            raise Exception(f"CoordinateD.__mul__({c1},{c2})")
        return x

    def __neg__(self):
        cdef CoordinateD x = CoordinateD()
        x.xx = x.xx - self.xx
        return x

    def __pos__(self):
        return self
