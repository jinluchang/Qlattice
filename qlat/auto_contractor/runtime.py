#    Qlattice (https://github.com/jinluchang/qlattice)
#
#    Copyright (C) 2023
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

from qlat import \
        timer, \
        timer_flops

from qlat_utils.ama import \
        ama_list, \
        ama_apply1, \
        ama_counts, \
        ama_extract

from qlat_utils import \
        load_prop, \
        get_gamma_matrix, \
        wilson_matrix_g5_herm, \
        mat_tr_sm, \
        mat_tr_cm, \
        mat_tr_wm, \
        mat_tr_wm_wm, \
        mat_tr_wm_sm, \
        mat_tr_sm_wm, \
        mat_tr_sm_sm, \
        mat_tr_wm_cm, \
        mat_tr_cm_wm, \
        mat_tr_cm_cm, \
        mat_mul_wm_wm, \
        mat_mul_wm_sm, \
        mat_mul_sm_wm, \
        mat_mul_sm_sm, \
        mat_mul_wm_cm, \
        mat_mul_cm_wm, \
        mat_mul_cm_cm

from . import auto_fac_funcs as aff
