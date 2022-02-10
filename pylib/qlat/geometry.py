import cqlat as c

class Geometry:

    def __init__(self, total_site, multiplicity = None):
        if multiplicity is None:
            self.cdata = c.mk_geo(total_site)
        else:
            self.cdata = c.mk_geo(total_site, multiplicity)

    def __del__(self):
        c.free_geo(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, Geometry)
        c.set_geo(self, v1)
        return self

    def copy(self, is_copying_data = True):
        x = Geometry([ 0, 0, 0, 0, ])
        if is_copying_data:
            x @= self
        return x

    def total_site(self):
        return c.get_total_site_geo(self)

    def total_volume(self):
        return c.get_total_volume_geo(self)

    def local_volume(self):
        return c.get_local_volume_geo(self)

    def multiplicity(self):
        return c.get_multiplicity_geo(self)

    def node_site(self):
        return c.get_node_site_geo(self)

    def eo(self):
        return c.get_eo_geo(self)

    def expansion_left(self):
        return c.get_expansion_left_geo(self)

    def expansion_right(self):
        return c.get_expansion_right_geo(self)

    def id_node(self):
        return c.get_id_node_geo(self)

    def num_node(self):
        return c.get_num_node_geo(self)

    def coor_node(self):
        return c.get_coor_node_geo(self)

    def size_node(self):
        return c.get_size_node_geo(self)

    def show(self):
        return "Geometry({}, {})".format(
                str(self.total_site()),
                self.multiplicity())

    def show_all(self):
        return "Geometry({}, {}, expansion_left={}, expansion_right={}, eo={})".format(
                str(self.total_site()),
                self.multiplicity(),
                str(self.expansion_left()),
                str(self.expansion_right()),
                self.eo())

    def coordinate_g_from_l(self, xl):
        return c.coordinate_g_from_l_geo(self, xl)

    def coordinate_l_from_g(self, xg):
        return c.coordinate_l_from_g_geo(self, xg)

    def is_local(self, xl):
        return c.is_local_geo(self, xl)

###

def geo_reform(geo, multiplicity = 1, expansion_left = None, expansion_right = None):
    assert isinstance(geo, Geometry)
    if expansion_left is None:
        expansion_left = [ 0, 0, 0, 0, ]
    elif isinstance(expansion_left, int):
        e = expansion_left
        expansion_left = [ e, e, e, e, ]
    if expansion_right is None:
        expansion_right = [ 0, 0, 0, 0, ]
    elif isinstance(expansion_right, int):
        e = expansion_right
        expansion_right = [ e, e, e, e, ]
    geo_new = Geometry([ 0, 0, 0, 0, ])
    c.set_geo_reform(geo_new, geo, multiplicity, expansion_left, expansion_right)
    return geo_new

def geo_eo(geo, eo = 0):
    assert isinstance(geo, Geometry)
    geo_new = Geometry([ 0, 0, 0, 0, ])
    c.set_geo_eo(geo_new, geo, eo)
    return geo_new
