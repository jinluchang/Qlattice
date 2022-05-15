import cqlat as c

class WilsonMatrix:

    # self.cdata

    def __init__(self):
        self.cdata = c.mk_wilson_matrix()

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_wilson_matrix(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, WilsonMatrix)
        c.set_wilson_matrix(self, v1)
        return self

    def copy(self, is_copying_data = True):
        x = WilsonMatrix()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __getstate__(self):
        # return a 12x12 2-D numpy array contain the data
        arr = np.frombuffer(c.get_state_wm(self), dtype = complex).reshape(12, 12)
        return arr

    def __setstate__(self, arr):
        c.set_state(self, arr.tobytes())

    def set_zero(self):
        return c.set_zero_wilson_matrix(self)

    def set_value(self, value_as_list_of_complex):
        return c.set_value_wilson_matrix(self, value_as_list_of_complex)

    def g5_herm(self):
        wm = WilsonMatrix()
        c.set_g5_herm_wilson_matrix(wm, self)
        return wm

###

class SpinMatrix:

    # self.cdata

    def __init__(self):
        self.cdata = c.mk_spin_matrix()

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_spin_matrix(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, SpinMatrix)
        c.set_spin_matrix(self, v1)
        return self

    def copy(self, is_copying_data = True):
        x = SpinMatrix()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_zero(self):
        return c.set_zero_spin_matrix(self)

    def set_value(self, value_as_list_of_complex):
        return c.set_value_spin_matrix(self, value_as_list_of_complex)

###
