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
            if other == 0:
                return 0
            return mat_mul_a_s(other, self)
        elif isinstance(other, SpinMatrix):
            return SpinMatrix(self.m @ other.m)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            if other == 0:
                return 0
            return mat_mul_a_s(other, self)
        elif isinstance(other, SpinMatrix):
            return SpinMatrix(other.m @ self.m)
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
            if other == 0:
                return 0
            return mat_mul_a_sc(other, self)
        elif isinstance(other, SpinMatrix):
            return mat_mul_sc_s(self, other)
        elif isinstance(other, SpinColorMatrix):
            return mat_mul_sc_sc(self, other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            if other == 0:
                return 0
            return mat_mul_a_sc(other, self)
        elif isinstance(other, SpinMatrix):
            return mat_mul_s_sc(other, self)
        elif isinstance(other, SpinColorMatrix):
            return mat_mul_sc_sc(other, self)
        else:
            return NotImplemented

    def trace(self):
        return self.m.trace()

###

def mat_mul_sc_sc(x, y):
    if x == 0 or y == 0:
        return 0
    return SpinColorMatrix(x.m @ y.m)

def mat_mul_sc_s(x, y):
    if x == 0:
        return 0
    return SpinColorMatrix(as_cont((as_cont(y.m.transpose()) @ x.m.reshape(12, 4, 3)).reshape(12, 12)))

def mat_mul_s_sc(x, y):
    if y == 0:
        return 0
    return SpinColorMatrix(as_cont((x.m @ y.m.reshape(4, 36)).reshape(12, 12)))

def mat_mul_s_s(x, y):
    return SpinMatrix(x.m @ y.m)

def mat_mul_a_s(coef, x):
    return SpinMatrix(x.m * coef)

def mat_mul_a_sc(coef, x):
    return SpinColorMatrix(x.m * coef)

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

def get_spin_matrix(op):
    assert op.otype == "G"
    assert op.s1 == "auto" and op.s2 == "auto"
    assert op.tag in [0, 1, 2, 3, 5]
    return gamma_matrix_list[op.tag]

def as_mspincolor(x):
    if isinstance(x, np.ndarray):
        return SpinColorMatrix(as_cont(x.reshape(12, 12)))
    elif x == 0:
        return SpinColorMatrix(np.zeros((12, 12,), dtype = complex))
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
