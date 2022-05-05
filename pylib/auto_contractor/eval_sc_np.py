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
import qlat as q

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
            return SpinMatrix(self.m * other)
        elif isinstance(other, SpinMatrix):
            return SpinMatrix(np.tensordot(self.m, other.m, 1))
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return SpinMatrix(other * self.m)
        elif isinstance(other, SpinMatrix):
            return SpinMatrix(np.tensordot(other.m, self.m, 1))
        else:
            return NotImplemented

    def trace(self):
        return np.trace(self.m)

###

class SpinColorMatrix:

    # self.m

    # self.m.shape = (4, 3, 3, 4,)

    def __init__(self, m):
        assert isinstance(m, np.ndarray)
        assert m.dtype == complex
        if m.shape == (4, 3, 3, 4,):
            self.m = m
        elif m.shape == (4, 4, 3, 3,):
            self.m = m.transpose(0, 2, 3, 1)
        elif m.shape == (144,):
            self.m = m.reshape(4, 4, 3, 3).transpose(0, 2, 3, 1)
        else:
            assert False

    def __mul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return SpinColorMatrix(self.m * other)
        elif isinstance(other, SpinMatrix):
            return SpinColorMatrix(np.tensordot(self.m, other.m, 1))
        elif isinstance(other, SpinColorMatrix):
            return SpinColorMatrix(np.tensordot(self.m, other.m, ((3, 2,), (0, 1,),)))
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return SpinMatrix(other * self.m)
        elif isinstance(other, SpinMatrix):
            return SpinColorMatrix(np.tensordot(other.m, self.m, 1))
        elif isinstance(other, SpinColorMatrix):
            return SpinColorMatrix(np.tensordot(other.m, self.m, ((3, 2,), (0, 1,),)))
        else:
            return NotImplemented

    def trace(self):
        return np.trace(self.m, axis1 = 0, axis2 = 3).trace()

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

def get_spin_matrix(op):
    assert op.otype == "G"
    assert op.s1 == "auto" and op.s2 == "auto"
    assert op.tag in [0, 1, 2, 3, 5]
    return gamma_matrix_list[op.tag]

def as_mspincolor(x):
    assert isinstance(x, np.ndarray)
    return SpinColorMatrix(x)

def adj_msc(x):
    return SpinColorMatrix(x.m.conj().transpose())

@q.timer
def g5_herm(x):
    x_h = gamma_matrix_list[5] * adj_msc(x) * gamma_matrix_list[5]
    return x_h

def msc_trace(x):
    return as_msc(x).trace()

def msc_trace2(x, y):
    if isinstance(x, SpinColorMatrix) and isinstance(y, SpinColorMatrix):
        v = np.tensordot(x.m, y.m, ((3, 2, 1, 0,), (0, 1, 2, 3,),)).item()
        return v
    elif isinstance(x, SpinColorMatrix) and isinstance(y, SpinMatrix):
        v = np.tensordot(x.m, y.m, ((3, 0,), (0, 1,),)).trace()
        return v
    elif isinstance(x, SpinMatrix) and isinstance(y, SpinColorMatrix):
        v = np.tensordot(x.m, y.m, ((1, 0,), (0, 3,),)).trace()
        return v
    return as_msc(x * y).trace()

def as_msc(x):
    assert isinstance(x, SpinColorMatrix)
    return x
