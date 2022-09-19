#    Qlattice (https://github.com/waterret/qlattice)
#
#    Copyright (C) 2021
#
#    Author: Luchang Jin (ljin.luchang@gmail.com)
#    Author: Masaaki Tomii
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import numpy as np
from numpy import ascontiguousarray as as_cont
import qlat as q
from auto_contractor.ama import *

class SpinMatrix:

    # self.m

    # self.m.shape = (4, 4,)

    def __init__(self, m):
        assert isinstance(m, np.ndarray)
        assert m.dtype == complex
        assert m.shape == (4, 4,)
        self.m = m

    def __mul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return SpinMatrix(mat_mul_a_s(other, self.m))
        elif isinstance(other, SpinMatrix):
            return SpinMatrix(mat_mul_s_s(self.m, other.m))
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return SpinMatrix(mat_mul_a_s(other, self.m))
        elif isinstance(other, SpinMatrix):
            return SpinMatrix(mat_mul_s_s(other.m, self.m))
        else:
            return NotImplemented

    def trace(self):
        return self.m.trace()

###

class SpinColorMatrix:

    # self.m

    # self.m.shape = (12, 12,)

    def __init__(self, m):
        assert isinstance(m, np.ndarray)
        assert m.dtype == complex
        assert m.shape == (12, 12,)
        self.m = m

    def __mul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return SpinColorMatrix(mat_mul_a_sc(other, self.m))
        elif isinstance(other, SpinMatrix):
            return SpinColorMatrix(mat_mul_sc_s(self.m, other.m))
        elif isinstance(other, SpinColorMatrix):
            return SpinColorMatrix(mat_mul_sc_sc(self.m, other.m))
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return SpinColorMatrix(mat_mul_a_sc(other, self.m))
        elif isinstance(other, SpinMatrix):
            return SpinColorMatrix(mat_mul_s_sc(other.m, self.m))
        elif isinstance(other, SpinColorMatrix):
            return SpinColorMatrix(mat_mul_sc_sc(other.m, self.m))
        else:
            return NotImplemented

    def trace(self):
        return self.m.trace()

###

def get_mat(x):
    return x.m

def mat_mul_sc_sc(x, y):
    return x @ y

def mat_mul_sc_s(x, y):
    # return as_cont((x.reshape(12, 4, 3).transpose(0, 2, 1) @ y).transpose(0, 2, 1).reshape(12, 12))
    return as_cont((as_cont(y.transpose()) @ x.reshape(12, 4, 3)).reshape(12, 12))

def mat_mul_s_sc(x, y):
    return as_cont((x @ y.reshape(4, 36)).reshape(12, 12))

def mat_mul_s_s(x, y):
    return x @ y

def mat_mul_a_s(coef, x):
    return x * coef

def mat_mul_a_sc(coef, x):
    return x * coef

def mat_sc_trace(x):
    return x.trace()

def mat_sc_sc_trace(x, y):
    return np.dot(x.ravel(), y.transpose().ravel()).item()

def mat_sc_s_trace(x, y):
    return np.tensordot(x.reshape(4, 3, 4, 3,), y, ((2, 0,), (0, 1,),)).trace()

def mat_s_sc_trace(x, y):
    return np.tensordot(x, y.reshape(4, 3, 4, 3,), ((1, 0,), (0, 2,),)).trace()

###

gamma_matrix_list = [
        # gamma_x
        SpinMatrix(
            (-1) * (-1j) *
            np.array([
                0, 0, 0, 1,
                0, 0, 1, 0,
                0, -1, 0, 0,
                -1, 0, 0, 0,
                ],
                dtype = complex).reshape(4, 4)
            ),
        # gamma_y
        SpinMatrix(
            (-1j) *
            np.array([
                0, 0, 0, -1j,
                0, 0, 1j, 0,
                0, 1j, 0, 0,
                -1j, 0, 0, 0,
                ],
                dtype = complex).reshape(4, 4)
            ),
        # gamma_z
        SpinMatrix(
            (-1) * (-1j) *
            np.array([
                0, 0, 1, 0,
                0, 0, 0, -1,
                -1, 0, 0, 0,
                0, 1, 0, 0,
                ],
                dtype = complex).reshape(4, 4)
            ),
        # gamma_t
        SpinMatrix(
            np.array([
                0, 0, 1, 0,
                0, 0, 0, 1,
                1, 0, 0, 0,
                0, 1, 0, 0,
                ],
                dtype = complex).reshape(4, 4)
            ),
        # gamma(4)
        None,
        # gamma_5 = gamma_x * gamma_y * gamma_z * gamma_t;
        SpinMatrix(
            np.array([
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, -1, 0,
                0, 0, 0, -1,
                ],
                dtype = complex).reshape(4, 4)
            ),
        ]

def get_gamma_matrix(mu):
    return gamma_matrix_list[mu]

def as_mspincolor(x):
    if isinstance(x, np.ndarray):
        return SpinColorMatrix(as_cont(x.reshape(12, 12)))
    elif x == 0:
        return SpinColorMatrix(np.zeros((12, 12,), dtype = complex))
    elif isinstance(x, q.WilsonMatrix):
        return SpinColorMatrix(as_cont(x.get_value()))
    elif isinstance(x, SpinColorMatrix):
        return x
    else:
        print("as_mspincolor:", x)
        assert False

def as_mspin(x):
    if isinstance(x, np.ndarray):
        return SpinMatrix(as_cont(x.reshape(4, 4)))
    elif x == 0:
        return SpinMatrix(np.zeros((4, 4,), dtype = complex))
    else:
        assert False

def adj_msc(x):
    return SpinColorMatrix(x.m.transpose().conj())

def g5_herm(x):
    if isinstance(x, (int, float)):
        return x
    elif isinstance(x, complex):
        return x.conjugate()
    else:
        x_h = adj_msc(x)
        vm = x_h.m.reshape(4, 3, 4, 3)
        assert vm.base is x_h.m
        vm[2:4, :, 0:2, :] *= -1
        vm[0:2, :, 2:4, :] *= -1
        return x_h

def msc_trace(x):
    assert isinstance(x, SpinColorMatrix)
    return x.trace()

def msc_trace2(x, y):
    if isinstance(x, SpinColorMatrix) and isinstance(y, SpinColorMatrix):
        return np.dot(x.m.ravel(), y.m.transpose().ravel()).item()
    elif isinstance(x, SpinColorMatrix) and isinstance(y, SpinMatrix):
        v = np.tensordot(x.m.reshape(4, 3, 4, 3,), y.m, ((2, 0,), (0, 1,),)).trace()
        return v
    elif isinstance(x, SpinMatrix) and isinstance(y, SpinColorMatrix):
        v = np.tensordot(x.m, y.m.reshape(4, 3, 4, 3,), ((1, 0,), (0, 2,),)).trace()
        return v
    else:
        assert False
        return as_msc(x * y).trace()

def as_msc(x):
    assert isinstance(x, SpinColorMatrix)
    return x
