#    Qlattice (https://github.com/jinluchang/qlattice)
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

import numpy as np

import functools

import qlat as q

@q.timer
def einsum_cache_path(subscripts, *operands, optimize = None):
    """
    Behave as ``np.einsum`` except will update ``optimize`` if it is an empty list
    """
    if optimize is None:
        return np.einsum(subscripts, *operands)
    if optimize == []:
        a, b = np.einsum_path(subscripts, *operands, optimize = 'optimal')
        optimize[:] = a
        q.displayln_info(a)
        q.displayln_info(b)
    return np.einsum(subscripts, *operands, optimize = optimize)

def as_wilson_matrix(x):
    if isinstance(x, np.ndarray):
        return x
    elif x == 0:
        return 0
    else:
        raise Exception(f"as_wilson_matrix: {x}")

einsum_optimize_g5_herm = []

def wilson_matrix_g5_herm(x):
    """
    corr[s1, n1, s2, n2]
    return gamma_5 * corr^dagger * gamma_5
    should equal to the original
    corr
    """
    g5 = get_gamma_matrix(5)
    corr_g5 = einsum_cache_path("ij,kbja,kl->ialb", g5, x.conj(), g5, optimize = einsum_optimize_g5_herm)
    return corr_g5

def as_wilson_matrix_g5_herm(x):
    if isinstance(x, np.ndarray):
        return wilson_matrix_g5_herm(x)
    elif x == 0:
        return 0
    else:
        raise Exception(f"as_wilson_matrix_g5_herm: {x}")

def load_prop(x):
    if isinstance(x, tuple):
        assert len(x) == 2 and x[0] == "g5_herm"
        return q.ama_apply1(as_wilson_matrix_g5_herm, x[1])
    return q.ama_apply1(as_wilson_matrix, x)

@functools.lru_cache(maxsize=None)
def get_gamma_matrix(mu):
    sm = q.get_gamma_matrix(mu)
    arr = np.asarray(sm)
    return arr

einsum_optimize_mat_tr_sm = []

def mat_tr_sm(mat):
    v = einsum_cache_path("ii->", mat, optimize = einsum_optimize_mat_tr_sm)
    return v

einsum_optimize_mat_tr_cm = []

def mat_tr_cm(mat):
    v = einsum_cache_path("ii->", mat, optimize = einsum_optimize_mat_tr_cm)
    return v

einsum_optimize_mat_tr_wm = []

def mat_tr_wm(mat):
    v = einsum_cache_path("iaia->", mat, optimize = einsum_optimize_mat_tr_wm)
    return v

einsum_optimize_mat_tr_wm_wm = []

def mat_tr_wm_wm(mat1, mat2):
    v = einsum_cache_path("iajb,jbia->", mat1, mat2, optimize = einsum_optimize_mat_tr_wm_wm)
    return v

einsum_optimize_mat_tr_wm_sm = []

def mat_tr_wm_sm(mat1, mat2):
    v = einsum_cache_path("iaja,ji->", mat1, mat2, optimize = einsum_optimize_mat_tr_wm_sm)
    return v

def mat_tr_sm_wm(mat1, mat2):
    return mat_tr_wm_sm(mat2, mat1)

einsum_optimize_mat_tr_sm_sm = []

def mat_tr_sm_sm(mat1, mat2):
    v = einsum_cache_path("ij,ji->", mat1, mat2, optimize = einsum_optimize_mat_tr_sm_sm)
    return v

einsum_optimize_mat_tr_wm_cm = []

def mat_tr_wm_cm(mat1, mat2):
    v = einsum_cache_path("iaib,ba->", mat1, mat2, optimize = einsum_optimize_mat_tr_wm_cm)
    return v

def mat_tr_cm_wm(mat1, mat2):
    return mat_tr_wm_cm(mat2, mat1)

einsum_optimize_mat_tr_cm_cm = []

def mat_tr_cm_cm(mat1, mat2):
    v = einsum_cache_path("ab,ba->", mat1, mat2, optimize = einsum_optimize_mat_tr_cm_cm)
    return v

einsum_optimize_mat_mul_wm_wm = []

def mat_mul_wm_wm(mat1, mat2):
    v = einsum_cache_path("iajb,jbkc->iakc", mat1, mat2, optimize = einsum_optimize_mat_mul_wm_wm)
    return v

einsum_optimize_mat_mul_wm_sm = []

def mat_mul_wm_sm(mat1, mat2):
    v = einsum_cache_path("iajb,jk->iakb", mat1, mat2, optimize = einsum_optimize_mat_mul_wm_sm)
    return v

einsum_optimize_mat_mul_sm_wm = []

def mat_mul_sm_wm(mat1, mat2):
    v = einsum_cache_path("ij,jakb->iakb", mat1, mat2, optimize = einsum_optimize_mat_mul_sm_wm)
    return v

einsum_optimize_mat_mul_sm_sm = []

def mat_mul_sm_sm(mat1, mat2):
    v = einsum_cache_path("ij,jk->ik", mat1, mat2, optimize = einsum_optimize_mat_mul_sm_sm)
    return v

einsum_optimize_mat_mul_wm_cm = []

def mat_mul_wm_cm(mat1, mat2):
    v = einsum_cache_path("iajb,bc->iajc", mat1, mat2, optimize = einsum_optimize_mat_mul_wm_cm)
    return v

einsum_optimize_mat_mul_cm_wm = []

def mat_mul_cm_wm(mat1, mat2):
    v = einsum_cache_path("ab,ibjc->iajc", mat1, mat2, optimize = einsum_optimize_mat_mul_cm_wm)
    return v

einsum_optimize_mat_mul_cm_cm = []

def mat_mul_cm_cm(mat1, mat2):
    v = einsum_cache_path("ab,bc->ac", mat1, mat2, optimize = einsum_optimize_mat_mul_cm_cm)
    return v
