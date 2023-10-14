from qlat_utils.complex cimport *
from qlat_utils.mat cimport *

cdef fused ElemType:
    ColorMatrix
    WilsonMatrix 
    NonRelWilsonMatrix
    IsospinMatrix
    SpinMatrix
    WilsonVector
    Complex
    ComplexF
    Double
    Float
    Long
    Int64t
    Int8t
    Char
