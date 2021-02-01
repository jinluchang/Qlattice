import cqlat as c

class Field:

    def __init__(self, ctype, geo = None, multiplicity = None):
        self.ctype = ctype
        if geo == None:
            self.cdata = c.mk_field(ctype)
        elif multiplicity == None:
            self.cdata = c.mk_field(ctype, geo.cdata)
        else:
            self.cdata = c.mk_field(ctype, geo.cdata, multiplicity)

    def __del__(self):
        c.free_field(self.ctype, self.cdata)

