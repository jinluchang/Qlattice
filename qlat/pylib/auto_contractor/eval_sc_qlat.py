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

import cqlat as c
import qlat as q
import numpy as np

from qlat import mat_mul_sc_sc, mat_mul_sc_s, mat_mul_s_sc, mat_mul_s_s, mat_mul_a_s, mat_mul_a_sc, mat_sc_trace, mat_sc_sc_trace, mat_sc_s_trace, mat_s_sc_trace

def as_mspin(x):
    if isinstance(x, q.SpinMatrix):
        return x
    sm = q.SpinMatrix()
    sm.set_value(x)
    return sm

def as_mspincolor(x):
    if isinstance(x, q.WilsonMatrix):
        return x
    wm = q.WilsonMatrix()
    if isinstance(x, np.ndarray):
        wm.set_value(x)
    elif x == 0:
        wm.set_zero()
    else:
        assert False
    return wm

def g5_herm(x):
    return x.g5_herm()

###

def get_mat(x):
    return x

gamma_matrix_list = [
        # gamma_x
        as_mspin(
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
        as_mspin(
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
        as_mspin(
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
        as_mspin(
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
        as_mspin(
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

###

def msc_trace(x):
    assert isinstance(x, q.WilsonMatrix)
    return mat_sc_trace(x, y)

def msc_trace2(x, y):
    if isinstance(x, q.WilsonMatrix) and isinstance(y, q.WilsonMatrix):
        return mat_sc_sc_trace(x, y)
    elif isinstance(x, q.WilsonMatrix) and isinstance(y, q.SpinMatrix):
        return mat_sc_s_trace(x, y)
    elif isinstance(x, q.SpinMatrix) and isinstance(y, q.WilsonMatrix):
        return mat_s_sc_trace(x, y)
    else:
        assert False
