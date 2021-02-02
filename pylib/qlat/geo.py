import cqlat as c

class Geometry:

    def __init__(self, total_site, multiplicity=None):
        if multiplicity == None:
            self.cdata = c.mk_geo(total_site)
        else:
            self.cdata = c.mk_geo(total_site, multiplicity)

    def __del__(self):
        c.free_geo(self.cdata)

    def total_site(self):
        return c.get_total_site_geo(self.cdata)

    def multiplicity(self):
        return c.get_multiplicity_geo(self.cdata)

    def node_site(self):
        return c.get_node_site_geo(self.cdata)

    def eo(self):
        return c.get_eo_geo(self.cdata)

    def expansion_left(self):
        return c.get_expansion_left_geo(self.cdata)

    def expansion_right(self):
        return c.get_expansion_right_geo(self.cdata)

    def id_node(self):
        return c.get_id_node_geo(self.cdata)

    def num_node(self):
        return c.get_num_node_geo(self.cdata)

    def coor_node(self):
        return c.get_coor_node_geo(self.cdata)

    def size_node(self):
        return c.get_size_node_geo(self.cdata)

def geo_reform(geo,
        multiplicity = 1,
        expansion_left = (0, 0, 0, 0),
        expansion_right = (0, 0, 0, 0)):
    if type(geo) == Geometry:
        geo = Geometry((0, 0, 0, 0))
        geo.cdata = c.geo_reform(geo.cdata, multiplicity, expansion_left, expansion_right)
        return geo
    else:
        raise Exception("geo_reform")

def geo_eo(geo, eo = 0):
    if type(geo) == Geometry:
        geo = Geometry((0, 0, 0, 0))
        geo.cdata = c.geo_eo(geo.cdata, eo)
        return geo
    else:
        raise Exception("geo_eo")

def show_geo(geo):
    if type(geo) == Geometry:
        return "Geometry({}, {})".format(
                str(geo.total_site()),
                geo.multiplicity())
    else:
        raise Exception("show_geo")

def show_geo_all(geo):
    if type(geo) == Geometry:
        return "Geometry({}, {}, expansion_left={}, expansion_right={}, eo={})".format(
                str(geo.total_site()),
                geo.multiplicity(),
                str(geo.expansion_left()),
                str(geo.expansion_right()),
                geo.eo())
    else:
        raise Exception("show_geo")

