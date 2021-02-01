import cqlat as c

class Geometry:

    def __init__(self, total_site, multiplicity=None):
        if multiplicity == None:
            self.cdata = c.mk_geo(total_site)
        else:
            self.cdata = c.mk_geo(total_site, multiplicity)

    def __del__(self):
        c.free_geo(self.cdata)


