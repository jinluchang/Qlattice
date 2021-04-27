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

    def set_zero(self):
        c.set_zero_lat_data(x)

    def set_dim_sizes(self, dim_sizes):
        c.set_dim_sizes_lat_data(self, dim_sizes)
