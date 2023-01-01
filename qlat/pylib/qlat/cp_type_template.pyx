### -------------------------------------------------------------------

cdef class FieldTYPENAME(FieldBase):

    ctype = ElemTypeTYPENAME

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __init__(self, geo = None, multiplicity = None):
        if geo is None:
            return
        assert isinstance(geo, Geometry)
        if multiplicity is None:
            self.xx.init((<Geometry>geo).xx)
        else:
            assert isinstance(multiplicity, int)
            self.xx.init((<Geometry>geo).xx, <int>multiplicity)

### -------------------------------------------------------------------

field_type_dict[ElemTypeTYPENAME] = FieldTYPENAME

### -------------------------------------------------------------------
