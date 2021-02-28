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

    def load(self, path):
        c.load_lat_data(self, path)

    def save(self, path):
        c.save_lat_data(self, path)

    def show(self):
        return c.show_lat_data(self)

    def set_zero(self):
        c.set_zero_lat_data(x)
