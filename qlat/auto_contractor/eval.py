#    Qlattice (https://github.com/jinluchang/qlattice)
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

from .compile import *
from qlat_utils.ama import *

from qlat_utils.c import \
        as_wilson_matrix, as_wilson_matrix_g5_herm

from . import auto_fac_funcs as aff

import numpy as np
import qlat as q
import copy
import cmath
import math
import importlib
import time
import os
import glob
import subprocess
import functools

class CCExpr:

    """
    self.cexpr_all
    self.module
    self.base_positions_dict
    self.cexpr_function_bare
    self.total_sloppy_flops
    self.expr_names
    self.positions
    """

    def __init__(self, cexpr_all, module, *, base_positions_dict=None):
        self.cexpr_all = cexpr_all
        self.module = module
        if base_positions_dict is None:
            base_positions_dict = {}
        self.base_positions_dict = base_positions_dict
        # module.cexpr_function(positions_dict, get_prop, is_ama_and_sloppy=False) => val as 1-D np.array
        self.cexpr_function_bare = module.cexpr_function
        self.total_sloppy_flops = module.total_sloppy_flops
        cexpr = self.cexpr_all["cexpr_optimized"]
        self.expr_names = [ name for name, expr in cexpr.named_exprs ]
        self.positions = cexpr.positions

    def get_expr_names(self):
        return self.expr_names

    def cexpr_function(self, positions_dict, get_prop, is_ama_and_sloppy=False):
        assert self.cexpr_function_bare is not None
        pd = self.base_positions_dict.copy()
        pd.update(positions_dict)
        return self.cexpr_function_bare(positions_dict=pd, get_prop=get_prop, is_ama_and_sloppy=is_ama_and_sloppy)

# -----

@q.timer
def cache_compiled_cexpr(
        calc_cexpr, path,
        *,
        is_cython=True,
        is_distillation=False,
        base_positions_dict=None,
        ):
    """
    Return an ``CCExpr`` created from ``cexpr = calc_cexpr()`` and cache the results.\n
    Save cexpr object in pickle format for future reuse.
    Generate python code and save for future reuse.
    Create CCExpr with loaded python/cython module.
    Return fully loaded ``ccexpr``.
    !!!Note that the module will not be reloaded if it has been loaded before!!!
    """
    fname = q.get_fname()
    if is_cython:
        path = path + "_cy"
    else:
        path = path + "_py"
    fn_pickle = path + "/cexpr_all.pickle"
    @q.timer
    def compile_cexpr_meson_setup():
        subprocess.run(["meson", "setup", "build"], cwd = path)
    @q.timer
    def compile_cexpr_meson_compile():
        subprocess.run(["meson", "compile", "-C", "build"], cwd = path)
        objs = glob.glob(f"{path}/build/cexpr_code.*.so")
        if len(objs) != 1:
            raise Exception(f"WARNING: compile_cexpr_meson_compile: {objs}")
    @q.timer
    def calc_compile_cexpr():
        q.timer_fork()
        def compile_cexpr():
            cexpr_original = calc_cexpr()
            content_original = display_cexpr(cexpr_original)
            q.qtouch_info(path + "/cexpr_original.txt", content_original)
            return cexpr_original
        cexpr_original = q.pickle_cache_call(
                compile_cexpr, path + "/cexpr_original.pickle", is_sync_node=False)
        def optimize():
            cexpr_optimized = copy.deepcopy(cexpr_original)
            cexpr_optimized.optimize()
            content_optimized = display_cexpr(cexpr_optimized)
            q.qtouch_info(path + "/cexpr_optimized.txt", content_optimized)
            return cexpr_optimized
        cexpr_optimized = q.pickle_cache_call(
                optimize, path + "/cexpr_optimized.pickle", is_sync_node=False)
        def gen_code():
            code_py = cexpr_code_gen_py(
                    cexpr_optimized,
                    is_cython = is_cython,
                    is_distillation = is_distillation)
            if is_cython:
                fn_py = path + "/cexpr_code.pyx"
            else:
                fn_py = path + "/cexpr_code.py"
            q.qtouch_info(fn_py, code_py)
            return code_py
        code_py = q.pickle_cache_call(
                gen_code, path + f"/cexpr_code.pickle", is_sync_node=False)
        if is_cython:
            q.qtouch_info(path + "/meson.build", meson_build_content)
            compile_cexpr_meson_setup()
            compile_cexpr_meson_compile()
        cexpr_all = dict()
        cexpr_all["cexpr_original"] = cexpr_original
        cexpr_all["cexpr_optimized"] = cexpr_optimized
        cexpr_all["code_py"] = code_py
        q.save_pickle_obj(cexpr_all, fn_pickle)
        q.timer_display()
        q.timer_merge()
        return cexpr_optimized
    if q.get_id_node() == 0 and not q.does_file_exist(fn_pickle):
        calc_compile_cexpr()
    q.sync_node()
    while not q.does_file_exist(fn_pickle):
        q.displayln(3, f"{fname}: Node {q.get_id_node()}: waiting for '{fn_pickle}'.")
        time.sleep(0.5)
    cexpr_all = q.load_pickle_obj(fn_pickle)
    q.displayln_info(1, f"{fname}: Loading '{path}'.")
    if is_cython:
        module = importlib.import_module((path + "/build/cexpr_code").replace("/", "."))
    else:
        module = importlib.import_module((path + "/cexpr_code").replace("/", "."))
    q.displayln_info(1, f"{fname}: Loaded '{path}'.")
    ccexpr = CCExpr(cexpr_all, module, base_positions_dict=base_positions_dict)
    return ccexpr

@q.timer
def get_expr_names(ccexpr : CCExpr):
    return ccexpr.get_expr_names()

@q.timer
def eval_cexpr(ccexpr : CCExpr, *, positions_dict, get_prop, is_ama_and_sloppy=False):
    """
    return 1 dimensional np.array
    cexpr can be cexpr object or can be a compiled function
    xg = positions_dict[position]
    mat_mspincolor = get_prop(flavor, xg_snk, xg_src)
    e.g. ("point-snk", [ 1, 2, 3, 4, ]) = positions_dict["x_1"]
    e.g. flavor = "l"
    e.g. xg_snk = ("point-snk", [ 1, 2, 3, 4, ])
    if is_ama_and_sloppy: return (val_ama, val_sloppy,)
    if not is_ama_and_sloppy: return val_ama
    Note:
    cexpr_function(positions_dict, get_prop, is_ama_and_sloppy=False) => val as 1-D np.array
    """
    return ccexpr.cexpr_function(positions_dict, get_prop, is_ama_and_sloppy)

@q.timer
def benchmark_eval_cexpr(
        cexpr : CCExpr,
        *,
        benchmark_size=10,
        benchmark_num=10,
        benchmark_num_ama=2,
        benchmark_rng_state=None,
        ):
    if benchmark_rng_state is None:
        benchmark_rng_state = q.RngState("benchmark_eval_cexpr")
    expr_names = get_expr_names(cexpr)
    n_expr = len(expr_names)
    # prop_dict = {}
    size = q.Coordinate([ 8, 8, 8, 16, ])
    positions_vars = []
    for pos in cexpr.positions:
        if pos == "size":
            continue
        if pos in aff.auto_fac_funcs_list:
            continue
        if pos in cexpr.base_positions_dict:
            continue
        positions_vars.append(pos)
    n_pos = len(positions_vars)
    positions = [
            ("point", benchmark_rng_state.split(f"positions {pos_idx}").c_rand_gen(size),)
            for pos_idx in range(n_pos)
            ]
    #
    def mk_pos_dict(k):
        positions_dict = {}
        positions_dict["size"] = size
        idx_list = q.random_permute(list(range(n_pos)), benchmark_rng_state.split(f"pos_dict {k}"))
        for pos, idx in zip(positions_vars, idx_list):
            positions_dict[pos] = positions[idx]
        return positions_dict
    positions_dict_list = [ mk_pos_dict(k) for k in range(benchmark_size) ]
    #
    @functools.lru_cache(maxsize=None)
    def mk_prop(flavor, pos_snk, pos_src):
        prop = make_rand_spin_color_matrix(benchmark_rng_state.split(f"prop {flavor} {pos_snk} {pos_src}"))
        prop_ama = make_rand_spin_color_matrix(benchmark_rng_state.split(f"prop ama {flavor} {pos_snk} {pos_src}"))
        ama_val = mk_ama_val(prop, pos_src, [ prop, prop_ama, ], [ 0, 1, ], [ 1.0, 0.5, ])
        return ama_val
    @functools.lru_cache(maxsize=None)
    def mk_prop_uu(tag, p, mu):
        uu = make_rand_color_matrix(benchmark_rng_state.split(f"prop U {tag} {p} {mu}"))
        return uu
    #
    def convert_pos(p):
        p_tag, p_val = p
        return p_tag, tuple(p_val.to_list())
    #
    @q.timer
    def get_prop(ptype, *args):
        if ptype == "U":
            tag, p, mu = args
            p = convert_pos(p)
            return mk_prop_uu(tag, p, mu)
        else:
            flavor = ptype
            pos_snk, pos_src = args
            pos_snk = convert_pos(pos_snk)
            pos_src = convert_pos(pos_src)
            return ama_extract(mk_prop(flavor, pos_snk, pos_src), is_sloppy=True)
    @q.timer
    def get_prop_ama(ptype, *args):
        if ptype == "U":
            tag, p, mu = args
            p = convert_pos(p)
            return mk_prop_uu(tag, p, mu)
        else:
            flavor = ptype
            pos_snk, pos_src = args
            pos_snk = convert_pos(pos_snk)
            pos_src = convert_pos(pos_src)
            return mk_prop(flavor, pos_snk, pos_src)
    #
    @q.timer_verbose
    def benchmark_eval_cexpr_run():
        res_list = []
        for k in range(benchmark_size):
            res = eval_cexpr(cexpr, positions_dict=positions_dict_list[k], get_prop=get_prop)
            res_list.append(res)
        res = np.array(res_list)
        assert res.shape == (benchmark_size, n_expr,)
        return res
    @q.timer_verbose
    def benchmark_eval_cexpr_run_with_ama():
        res_list = []
        for k in range(benchmark_size):
            res1 = eval_cexpr(cexpr, positions_dict=positions_dict_list[k], get_prop=get_prop_ama)
            res2 = eval_cexpr(cexpr, positions_dict=positions_dict_list[k], get_prop=get_prop)
            res_ama, res_sloppy = eval_cexpr(cexpr, positions_dict=positions_dict_list[k], get_prop=get_prop_ama, is_ama_and_sloppy=True)
            assert q.qnorm(res1 - res_ama) == 0
            assert q.qnorm(res2 - res_sloppy) == 0
            res_list.append(res_ama)
        res = np.array(res_list)
        assert res.shape == (benchmark_size, n_expr,)
        return res
    def mk_check_vector(k):
        rs = benchmark_rng_state.split(f"check_vector {k}")
        res = np.array([
            [ complex(rs.u_rand_gen(1.0, -1.0), rs.u_rand_gen(1.0, -1.0)) for i in range(n_expr) ]
            for k in range(benchmark_size) ])
        return res
    check_vector_list = [ mk_check_vector(k) for k in range(3) ]
    def check_res(res):
        return [ np.tensordot(res, cv).item() for cv in check_vector_list ]
    q.displayln_info(f"benchmark_eval_cexpr: benchmark_size={benchmark_size}")
    q.timer_fork(0)
    check = None
    for i in range(benchmark_num):
        res = benchmark_eval_cexpr_run()
        new_check = check_res(res)
        if check is None:
            check = new_check
        else:
            assert check == new_check
    check_ama = None
    for i in range(benchmark_num_ama):
        res_ama = benchmark_eval_cexpr_run_with_ama()
        new_check_ama = check_res(res_ama)
        if check_ama is None:
            check_ama = new_check_ama
        else:
            assert check_ama == new_check_ama
    q.timer_display()
    q.timer_merge()
    q.displayln_info(f"benchmark_eval_cexpr: {benchmark_show_check(check)} {benchmark_show_check(check_ama)}")
    return check, check_ama

def get_mpi_chunk(total_list, *, rng_state = None):
    """
    rng_state has to be the same on all the nodes
    e.g. rng_state = q.RngState("get_mpi_chunk")
    """
    chunk_number = q.get_num_node()
    chunk_id = q.get_id_node()
    chunk_list = q.get_chunk_list(total_list, chunk_number = chunk_number, rng_state = rng_state)
    if chunk_id < len(chunk_list):
        return chunk_list[chunk_id]
    else:
        return []

# -----------------------------------------

def get_cexpr_names(ccexpr : CCExpr):
    q.displayln_info("WARNING: get_expr_names: use get_expr_names instead.")
    return get_expr_names(ccexpr)

meson_build_content = r"""project(
  'qlat-auto-contractor-cexpr', 'cpp', 'cython',
  version: '1.0',
  license: 'GPL-3.0-or-later',
  default_options: [
    'warning_level=3',
    'cpp_std=c++14',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    'cython_language=cpp',
    ])
#
add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])
#
qlat_utils_cpp = meson.get_compiler('cpp')
#
qlat_utils_py3 = import('python').find_installation('python3')
message(qlat_utils_py3.path())
message(qlat_utils_py3.get_install_dir())
#
qlat_utils_omp = dependency('openmp').as_system()
qlat_utils_zlib = dependency('zlib').as_system()
#
qlat_utils_math = qlat_utils_cpp.find_library('m')
#
qlat_utils_numpy_include = run_command(qlat_utils_py3, '-c', 'import numpy as np ; print(np.get_include())',
  check: true).stdout().strip()
message('numpy include', qlat_utils_numpy_include)
#
qlat_utils_numpy = declare_dependency(
  include_directories:  include_directories(qlat_utils_numpy_include),
  dependencies: [ qlat_utils_py3.dependency(), ],
  ).as_system()
#
if qlat_utils_cpp.check_header('Eigen/Eigen')
  qlat_utils_eigen = dependency('', required: false)
elif qlat_utils_cpp.check_header('Grid/Eigen/Eigen')
  qlat_utils_eigen = dependency('', required: false)
else
  qlat_utils_eigen = dependency('eigen3').as_system()
endif
#
qlat_utils_include = run_command(qlat_utils_py3, '-c', 'import qlat_utils as q ; print("\\n".join(q.get_include_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat_utils include', qlat_utils_include)
#
qlat_utils_lib = run_command(qlat_utils_py3, '-c', 'import qlat_utils as q ; print("\\n".join(q.get_lib_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat_utils lib', qlat_utils_lib)
#
qlat_utils_pxd = run_command(qlat_utils_py3, '-c', 'import qlat_utils as q ; print("\\n".join(q.get_pxd_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
# message('qlat_utils pxd', qlat_utils_pxd)
qlat_utils_pxd = files(qlat_utils_pxd)
#
qlat_utils_header = run_command(qlat_utils_py3, '-c', 'import qlat_utils as q ; print("\\n".join(q.get_header_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
# message('qlat_utils header', qlat_utils_header)
qlat_utils_header = files(qlat_utils_header)
#
qlat_utils = declare_dependency(
  include_directories: include_directories(qlat_utils_include),
  dependencies: [
    qlat_utils_py3.dependency().as_system(),
    qlat_utils_cpp.find_library('qlat-utils', dirs: qlat_utils_lib),
    qlat_utils_numpy, qlat_utils_eigen, qlat_utils_omp, qlat_utils_zlib, qlat_utils_math, ],
  )
#
py3 = import('python').find_installation('python3', pure: false)
#
deps = [ qlat_utils, ]
incdir = []
#
codelib = py3.extension_module('cexpr_code',
  files('cexpr_code.pyx'),
  dependencies: deps,
  include_directories: incdir,
  install: false,
  )
"""

def make_rand_spin_color_matrix(rng_state):
    rs = rng_state
    wm = q.WilsonMatrix()
    wm_arr = np.asarray(wm)
    wm_arr[:] = np.array(
            [ rs.u_rand_gen() + 1j * rs.u_rand_gen() for i in range(144) ],
            dtype = complex).reshape(12, 12)
    return wm

def make_rand_spin_matrix(rng_state):
    rs = rng_state
    sm = q.SpinMatrix()
    sm_arr = np.asarray(sm)
    sm_arr[:] = np.array(
            [ rs.u_rand_gen() + 1j * rs.u_rand_gen() for i in range(16) ],
            dtype = complex).reshape(4, 4)
    return sm

def make_rand_color_matrix(rng_state):
    rs = rng_state
    cm = q.ColorMatrix()
    cm_arr = np.asarray(cm)
    cm_arr[:] = np.array(
            [ rs.u_rand_gen() + 1j * rs.u_rand_gen() for i in range(9) ],
            dtype = complex).reshape(3, 3)
    return cm

def benchmark_show_check(check):
    return " ".join([ f"{v:.10E}" for v in check ])

def sqr_component(x):
    return x.real * x.real + 1j * x.imag * x.imag

def sqrt_component(x):
    return math.sqrt(x.real) + 1j * math.sqrt(x.imag)

def sqr_component_array(arr):
    return np.array([ sqr_component(x) for x in arr ])

def sqrt_component_array(arr):
    return np.array([ sqrt_component(x) for x in arr ])
