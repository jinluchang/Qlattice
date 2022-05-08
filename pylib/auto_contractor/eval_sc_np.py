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
            return SpinMatrix(self.m * other)
        elif isinstance(other, SpinMatrix):
            return SpinMatrix(self.m @ other.m)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            if other == 0:
                return 0
            return SpinMatrix(other * self.m)
        elif isinstance(other, SpinMatrix):
            return SpinMatrix(other.m @ self.m)
        else:
            return NotImplemented

    def trace(self):
        return self.m.trace()

###

class SpinColorMatrix:

    # self.m

    # self.m.shape = (4, 3, 3, 4,)

    def __init__(self, m):
        assert isinstance(m, np.ndarray)
        assert m.dtype == complex
        if m.shape == (12, 12,):
            self.m = m
        elif m.shape == (144,):
            self.m = m.reshape(12, 12)
        else:
            assert False

    def __mul__(self, other):
        if isinstance(other, (int, float, complex,)):
            if other == 0:
                return 0
            return SpinColorMatrix(self.m * other)
        elif isinstance(other, SpinMatrix):
            return SpinColorMatrix((other.m.transpose() @ self.m.reshape(12, 4, 3)).reshape(12, 12).copy())
        elif isinstance(other, SpinColorMatrix):
            return SpinColorMatrix(self.m @ other.m)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            if other == 0:
                return 0
            return SpinMatrix(other * self.m)
        elif isinstance(other, SpinMatrix):
            return SpinColorMatrix((other.m @ self.m.reshape(4, 36)).reshape(12, 12).copy())
        elif isinstance(other, SpinColorMatrix):
            return SpinColorMatrix(other.m @ self.m)
        else:
            return NotImplemented

    def trace(self):
        return self.m.trace()

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
        return SpinColorMatrix(x)
    elif isinstance(x, (int, float, complex)):
        return x
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
    if isinstance(x, SpinColorMatrix):
        return x.trace()
    elif isinstance(x, (int, float, complex)):
        if x == 0:
            return 0
        else:
            return x * 12
    else:
        assert False
        return as_msc(x).trace()

def msc_trace2(x, y):
    if isinstance(x, SpinColorMatrix) and isinstance(y, SpinColorMatrix):
        v = np.tensordot(x.m, y.m, ((1, 0,), (0, 1,),)).item()
        return v
    elif isinstance(x, SpinColorMatrix) and isinstance(y, SpinMatrix):
        v = np.tensordot(x.m.reshape(4, 3, 4, 3,), y.m, ((2, 0,), (0, 1,),)).trace()
        return v
    elif isinstance(x, SpinMatrix) and isinstance(y, SpinColorMatrix):
        v = np.tensordot(x.m, y.m.reshape(4, 3, 4, 3,), ((1, 0,), (0, 2,),)).trace()
        return v
    elif isinstance(x, (int, float, complex)):
        return x * msc_trace(y)
    elif isinstance(y, (int, float, complex)):
        return y * msc_trace(x)
    else:
        assert False
        return as_msc(x * y).trace()

def as_msc(x):
    assert isinstance(x, SpinColorMatrix)
    return x
