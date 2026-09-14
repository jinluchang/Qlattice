# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.fermion_action``
==============================\n
Fermion action parameterisation for Mobius and ZMobius domain-wall fermions.\n
Documentation: ``docs/qlat/qlat_fermion_action.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc

from cpython.long cimport PyLong_FromVoidPtr
from cpython.long cimport PyLong_AsVoidPtr

import cqlat as c

cdef inline cc.FermionAction* get_fermion_action_ptr(object fa) except? NULL:
    return <cc.FermionAction*>PyLong_AsVoidPtr(fa.cdata)

cdef inline cc.ComplexD to_complex_d(object z):
    cdef cc.PyComplexD pz = <cc.PyComplexD>z
    return (<cc.ComplexD*>(&pz))[0]

def mk_fermion_action_mobius(
        cc.RealD mass, cc.Int ls, cc.RealD m5, cc.RealD mobius_scale):
    cdef cc.FermionAction* pfa = new cc.FermionAction(
        mass, ls, m5, mobius_scale, True, False)
    return PyLong_FromVoidPtr(<void*>pfa)

def mk_fermion_action_zmobius(cc.RealD mass, cc.RealD m5, omega):
    cdef cc.Int ls = <cc.Int>len(omega)
    cdef cc.FermionAction* pfa = new cc.FermionAction(
        mass, ls, m5, 0.0, True, True)
    cdef Py_ssize_t i
    cdef object z
    cdef object b
    for i in range(ls):
        z = complex(omega[i])
        b = 0.5 * (1.0 / z + 1.0)
        pfa.bs[i] = to_complex_d(b)
        pfa.cs[i] = to_complex_d(b - 1.0)
    return PyLong_FromVoidPtr(<void*>pfa)

def get_mass_fermion_action(fa):
    cdef cc.FermionAction* pfa = get_fermion_action_ptr(fa)
    return pfa.mass

def get_m5_fermion_action(fa):
    cdef cc.FermionAction* pfa = get_fermion_action_ptr(fa)
    return pfa.m5

class FermionAction:
    def __init__(self, *, mass, ls, m5, mobius_scale=1.0, omega=None):
        assert isinstance(mass, float)
        assert isinstance(ls, int)
        assert isinstance(m5, float)
        if omega is None:
            self.cdata = mk_fermion_action_mobius(mass, ls, m5, mobius_scale)
        else:
            assert isinstance(omega, list)
            assert ls == len(omega)
            self.cdata = mk_fermion_action_zmobius(mass, m5, omega)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_fermion_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, FermionAction)
        c.set_fermion_action(self, v1)
        return self

    def mass(self):
        return get_mass_fermion_action(self)

    def ls(self):
        return c.get_ls_fermion_action(self)

    def m5(self):
        return get_m5_fermion_action(self)

    def omega(self):
        return c.get_omega_fermion_action(self)

    def mobius_scale(self):
        return c.get_mobius_scale_fermion_action(self)
