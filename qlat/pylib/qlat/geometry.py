import qlat.cqlat as c

class Geometry:

    # self.cdata

    def __init__(self, total_site = None, multiplicity = None):
        # if total_site is None: create geo uninitialized
        # elif multiplicity is None: create geo with multiplicity = 1
        self.cdata = c.mk_geo()
        if total_site is not None:
            if multiplicity is None:
                c.set_geo_total_site(self, total_site)
            else:
                c.set_geo_total_site(self, total_site, multiplicity)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_geo(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, Geometry)
        c.set_geo(self, v1)
        return self

    def copy(self, is_copying_data = True):
        x = Geometry()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

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

    def is_local_xg(self, xg):
        # return a global coordinate is inside the local volume or not
        return c.is_local_xg_geo(self, xg)

    def xg_list(self):
        # return xg for all local sites
        return c.get_xg_list(self)

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
    geo_new = geo.copy()
    c.set_geo_reform(geo_new, multiplicity, expansion_left, expansion_right)
    return geo_new

def geo_eo(geo, eo = 0):
    assert isinstance(geo, Geometry)
    geo_new = Geometry([ 0, 0, 0, 0, ])
    geo_new = geo.copy()
    c.set_geo_eo(geo_new, eo)
    return geo_new
