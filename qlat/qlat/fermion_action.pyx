# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.fermion_action``
==============================\n
Fermion action parameterisation for Mobius and ZMobius domain-wall fermions.
``FermionAction`` is a ``cdef class`` that owns its C++ object by value.\n
Documentation: ``docs/qlat/qlat_fermion_action.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc

cdef inline cc.ComplexD to_complex_d(object z):
    cdef cc.PyComplexD pz = <cc.PyComplexD>z
    return (<cc.ComplexD*>(&pz))[0]

cdef class FermionAction:

    def __cinit__(self):
        self.xx = cc.FermionAction()

    def __init__(self, *, mass, ls, m5, mobius_scale=1.0, omega=None):
        cdef Py_ssize_t i
        cdef object z
        cdef object b
        assert isinstance(mass, float)
        assert isinstance(ls, int)
        assert isinstance(m5, float)
        if omega is None:
            self.xx = cc.FermionAction(mass, ls, m5, mobius_scale, True, False)
        else:
            assert isinstance(omega, list)
            assert ls == len(omega)
            self.xx = cc.FermionAction(mass, ls, m5, 0.0, True, True)
            for i in range(ls):
                z = complex(omega[i])
                b = 0.5 * (1.0 / z + 1.0)
                self.xx.bs[i] = to_complex_d(b)
                self.xx.cs[i] = to_complex_d(b - 1.0)

    def __imatmul__(self, FermionAction v1):
        self.xx = v1.xx
        return self

    def mass(self):
        return self.xx.mass

    def ls(self):
        return self.xx.ls

    def m5(self):
        return self.xx.m5

    def omega(self):
        return get_omega_fermion_action(self)

    def mobius_scale(self):
        return self.xx.mobius_scale

### -------------------------------------------------------------------
### cqlat-compatible entry points

def get_mass_fermion_action(FermionAction fa):
    return fa.xx.mass

def get_m5_fermion_action(FermionAction fa):
    return fa.xx.m5

def set_fermion_action(FermionAction fa_new, FermionAction fa):
    fa_new.xx = fa.xx

def get_ls_fermion_action(FermionAction fa):
    return fa.xx.ls

def get_omega_fermion_action(FermionAction fa):
    if not fa.xx.is_using_zmobius:
        return None
    cdef cc.Long i
    cdef cc.std_vector[cc.RealD] re = cc.py_get_omega_fermion_action_re(fa.xx)
    cdef cc.std_vector[cc.RealD] im = cc.py_get_omega_fermion_action_im(fa.xx)
    cdef list omega = []
    for i in range(re.size()):
        omega.append(complex(re[i], im[i]))
    return omega

def get_mobius_scale_fermion_action(FermionAction fa):
    if fa.xx.is_using_zmobius:
        assert fa.xx.mobius_scale == 0.0
    else:
        assert fa.xx.mobius_scale != 0.0
    return fa.xx.mobius_scale
