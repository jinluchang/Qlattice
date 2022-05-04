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

import gpt as g
import numpy as np

def get_spin_matrix(op):
    assert op.otype == "G"
    assert op.s1 == "auto" and op.s2 == "auto"
    assert op.tag in [0, 1, 2, 3, 5]
    return g.gamma[op.tag]

def as_msc(x):
    return ascontiguoustensor(x)

def ascontiguoustensor(x):
    # isinstance(x, g.core.tensor)
    return g.tensor(np.ascontiguousarray(x.array), x.otype)

def as_mspincolor(x):
    if isinstance(x, g.core.tensor):
        return ascontiguoustensor(x)
    else:
        return g.tensor(np.ascontiguousarray(np.array(x)), g.ot_matrix_spin_color(4, 3))

def adj_msc(x):
    # isinstance(x, g.core.tensor)
    x = g.adj(x)
    return ascontiguoustensor(x)

def g5_herm(x):
    x_h = ascontiguoustensor(
            ascontiguoustensor(
                g.gamma[5]
                * adj_msc(x))
            * g.gamma[5])
    return x_h

def msc_trace(x):
    return g.trace(x)

def msc_trace2(x, y):
    return msc_trace(x * y)
