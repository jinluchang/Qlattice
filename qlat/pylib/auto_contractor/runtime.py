#    Qlattice (https://github.com/waterret/qlattice)
#
#    Copyright (C) 2022
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

from auto_contractor.eval import \
        load_prop, \
        ama_apply, ama_apply1, ama_apply2, ama_apply2_r, ama_apply2_l, ama_list, ama_counts, ama_extract

from qlat_utils import \
        rel_mod, rel_mod_sym, c_rel_mod_sqr

from qlat_utils.c import \
        WilsonMatrix, \
        SpinMatrix, \
        ColorMatrix, \
        get_gamma_matrix, \
        mat_tr_sm, \
        mat_tr_wm, \
        mat_tr_wm_wm, \
        mat_tr_wm_sm, \
        mat_tr_sm_wm, \
        mat_tr_sm_sm, \
        mat_tr_cm_wm, \
        mat_tr_wm_cm, \
        mat_tr_cm_cm, \
        mat_mul_wm_wm, \
        mat_mul_wm_sm, \
        mat_mul_sm_wm, \
        mat_mul_sm_sm, \
        mat_mul_wm_cm, \
        mat_mul_cm_wm, \
        mat_mul_cm_cm

from qlat import \
        timer, timer_flops

import numpy as np
