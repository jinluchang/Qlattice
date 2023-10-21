from . cimport everything as cqlat
from .field_types cimport FieldWilsonMatrix, FieldWilsonVector
from .selected_field_types cimport SelectedFieldWilsonMatrix
from .selected_points_types cimport SelectedPointsWilsonMatrix

cdef class Prop(FieldWilsonMatrix):

    cdef cqlat.Handle[cqlat.Prop] xxx(self)

cdef class SelProp(SelectedFieldWilsonMatrix):

    cdef cqlat.Handle[cqlat.SelProp] xxx(self)

cdef class PselProp(SelectedPointsWilsonMatrix):

    cdef cqlat.Handle[cqlat.PselProp] xxx(self)

cdef class FermionField4d(FieldWilsonVector):

    cdef cqlat.Handle[cqlat.FermionField4d] xxx(self)
