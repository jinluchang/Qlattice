from . cimport everything as cc

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
        return f"Coordinate({self.to_list()})"

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
        self.xx = cc.Coordinate(x[0], x[1], x[2], x[3])

    def __getstate__(self):
        return self.to_list()

    def __setstate__(self, state):
        self.from_list(state)

    def sqr(self):
        """
        Return the square sum of all the components as ``long``.
        """
        return cc.sqr(self.xx)

    def r_sqr(self):
        """
        get spatial distance square as int
        """
        return self.xx[0] * self.xx[0] + self.xx[1] * self.xx[1] + self.xx[2] * self.xx[2]

    def __getitem__(self, int key):
        assert 0 <= key
        assert key < 4
        return self.xx[key]

    def __setitem__(self, int key, int val):
        assert 0 <= key
        assert key < 4
        cdef int* p_val = &self.xx[key]
        p_val[0] = val

    def __add__(Coordinate c1, Coordinate c2):
        cdef Coordinate x = Coordinate()
        x.xx = c1.xx + c2.xx
        return x

    def __sub__(Coordinate c1, Coordinate c2):
        cdef Coordinate x = Coordinate()
        x.xx = c1.xx - c2.xx
        return x

    def __mul__(c1, c2):
        cdef Coordinate x = Coordinate()
        if isinstance(c2, int):
            x.xx = (<Coordinate>c1).xx * (<int>c2)
        elif isinstance(c2, Coordinate):
            x.xx = (<Coordinate>c1).xx * (<Coordinate>c2).xx
        else:
            raise Exception(f"Coordinate.__mul__({c1},{c2})")
        return x

    def __rmul__(c1, c2):
        return c1 * c2

    def __neg__(self):
        cdef Coordinate x = Coordinate()
        x.xx = x.xx - self.xx
        return x

    def __pos__(self):
        return self

    def __eq__(self, other):
        return isinstance(other, Coordinate) and self.xx == (<Coordinate>other).xx

    def from_index(self, long index, Coordinate size):
        self.xx = cc.coordinate_from_index(index, size.xx)

    def to_index(self, Coordinate size):
        return cc.index_from_coordinate(self.xx, size.xx)

# ------

def mod(Coordinate c, Coordinate size):
    """
    mod based on ``size``
    return ``x``
    ``0 <= x < size``
    """
    cdef Coordinate x = Coordinate()
    x.xx = cc.mod(c.xx, size.xx)
    return x

def smod(Coordinate c, Coordinate size):
    """
    smod based on ``size``
    return ``x``
    ``-size/2 <= x < size/2``
    """
    cdef Coordinate x = Coordinate()
    x.xx = cc.smod(c.xx, size.xx)
    return x

def middle_mod(Coordinate x, Coordinate y, Coordinate size):
    """
    return middle of x and y
    xm = mod(x, size)
    ym = mod(y, size)
    if xm<=ym: return mod(xm + smod(ym - xm, size), size)
    else: return mod(ym + smod(xm - ym, size), size)
    """
    cdef Coordinate ret = Coordinate()
    ret.xx = cc.middle_mod(x.xx, y.xx, size.xx)
    return ret

def coordinate_from_index(long index, size):
    cdef Coordinate x = Coordinate()
    if isinstance(size, Coordinate):
        x.xx = cc.coordinate_from_index(index, (<Coordinate>size).xx)
        return x
    else:
        x.xx = cc.coordinate_from_index(index, cc.Coordinate(size[0], size[1], size[2], size[3]))
        return x.to_list()

def index_from_coordinate(x, size):
    if isinstance(x, Coordinate) and isinstance(size, Coordinate):
        return cc.index_from_coordinate((<Coordinate>x).xx, (<Coordinate>size).xx)
    else:
        return cc.index_from_coordinate(cc.Coordinate(x[0], x[1], x[2], x[3]),
                                        cc.Coordinate(size[0], size[1], size[2], size[3]))

# ------

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
        if isinstance(c2, (int, float,)):
            x.xx = (<CoordinateD>c1).xx * (<double>c2)
        elif isinstance(c2, CoordinateD):
            x.xx = (<CoordinateD>c1).xx * (<CoordinateD>c2).xx
        else:
            raise Exception(f"CoordinateD.__mul__({c1},{c2})")
        return x

    def __rmul__(c1, c2):
        return c1 * c2

    def __neg__(self):
        cdef CoordinateD x = CoordinateD()
        x.xx = x.xx - self.xx
        return x

    def __pos__(self):
        return self

    def __eq__(self, other):
        return isinstance(other, CoordinateD) and self.xx == (<CoordinateD>other).xx

### ------------------------------------------------
