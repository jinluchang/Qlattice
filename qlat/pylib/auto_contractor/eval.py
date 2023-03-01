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

from auto_contractor.compile import *
from auto_contractor.ama import *

from qlat_utils import \
        rel_mod, rel_mod_sym, c_rel_mod_sqr

from qlat_utils.cp import \
        as_wilson_matrix, as_wilson_matrix_g5_herm

import numpy as np
import qlat as q
import copy
import cmath
import math
import importlib
import time
import os
import subprocess

def load_prop(x):
    if isinstance(x, tuple):
        assert len(x) == 2 and x[0] == "g5_herm"
        return ama_apply1(as_wilson_matrix_g5_herm, x[1])
    return ama_apply1(as_wilson_matrix, x)

@q.timer
def get_cexpr_names(cexpr):
    names = [ name for name, expr in cexpr.named_exprs ]
    return names

def eval_cexpr_get_props(cexpr : CExpr, *, positions_dict, get_prop):
    assert cexpr.function is not None
    return cexpr.function["cexpr_function_get_prop"](positions_dict, get_prop)

def eval_cexpr_eval(cexpr : CExpr, *, positions_dict, props):
    assert cexpr.function is not None
    return cexpr.function["cexpr_function_eval"](positions_dict, props)

def eval_cexpr(cexpr : CExpr, *, positions_dict, get_prop, is_ama_and_sloppy = False):
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
    """
    assert cexpr.function is not None
    return cexpr.function["cexpr_function"](positions_dict = positions_dict, get_prop = get_prop, is_ama_and_sloppy = is_ama_and_sloppy)

meson_build_content = """\
project('qlat-auto-contractor-cexpr', 'cpp', 'cython',
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
add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])
cpp = meson.get_compiler('cpp')
py_mod = import('python')
py3 = py_mod.find_installation('python3')
message(py3.path())
message(py3.get_install_dir())
omp = dependency('openmp').as_system()
zlib = dependency('zlib').as_system()
fftw = dependency('fftw3').as_system()
fftwf = dependency('fftw3f').as_system()
message('fftw libdir', fftw.get_variable('libdir'))
message('fftwf libdir', fftwf.get_variable('libdir'))
fftw_all = [ fftw, fftwf, ]
math = cpp.find_library('m')
numpy_include = run_command(py3.path(), '-c', 'import numpy as np ; print(np.get_include())', check: true).stdout().strip()
message('numpy include', numpy_include)
numpy = declare_dependency(
  include_directories:  include_directories(numpy_include),
  dependencies: [ py3.dependency(), ],
  ).as_system()
qlat_utils_include = run_command(py3.path(), '-c', 'import qlat_utils as q ; print("\\\\n".join(q.get_include_list()))', env: environment({'q_verbose': '-1'}), check: true).stdout().strip().split('\\n')
message('qlat_utils include', qlat_utils_include)
qlat_utils_lib = run_command(py3.path(), '-c', 'import qlat_utils as q ; print("\\\\n".join(q.get_lib_list()))', env: environment({'q_verbose': '-1'}), check: true).stdout().strip().split('\\n')
message('qlat_utils lib', qlat_utils_lib)
qlat_utils = declare_dependency(
  include_directories:  include_directories(qlat_utils_include),
  dependencies: [
    py3.dependency().as_system(),
    cpp.find_library('qlat-utils', dirs: qlat_utils_lib),
    numpy, omp, zlib, math, ],
  )
deps = [ qlat_utils, fftw_all, ]
if not cpp.check_header('Eigen/Eigen')
  eigen = dependency('eigen3').as_system()
  deps += [ eigen, ]
endif
incdir = []
codelib = py3.extension_module('code',
  files('code.pyx'),
  dependencies: deps,
  include_directories: incdir,
  install: false,
  )
"""

@q.timer
def cache_compiled_cexpr(calc_cexpr, path):
    """
    Return fully loaded ``cexpr = calc_cexpr()`` and cache the results\n
    Save cexpr object in pickle format for future reuse.
    Generate python code and save for future reuse.
    Load the python module and assign ``cexpr.function`` and ``cexpr.total_sloppy_flops``.
    Return fully loaded ``cexpr``.
    """
    fn_pickle = path + "/data.pickle"
    fn_code = path + "/build/code.so"
    @q.timer
    def compile_cexpr_meson_setup():
        subprocess.run(["meson", "setup", "build"], cwd = path)
    @q.timer
    def compile_cexpr_meson_compile():
        subprocess.run(["meson", "compile", "-C", "build"], cwd = path)
    @q.timer
    def calc_compile_cexpr():
        cexpr_original = calc_cexpr()
        cexpr_optimized = copy.deepcopy(cexpr_original)
        cexpr_optimized.optimize()
        code_py = cexpr_code_gen_py(cexpr_optimized)
        fn_py = path + "/code.pyx"
        q.qtouch_info(fn_py, code_py)
        content_original = display_cexpr(cexpr_original)
        q.qtouch_info(path + "/cexpr.original.txt", content_original)
        content_optimized = display_cexpr(cexpr_optimized)
        q.qtouch_info(path + "/cexpr.optimized.txt", content_optimized)
        q.qtouch_info(path + "/meson.build", meson_build_content)
        compile_cexpr_meson_setup()
        compile_cexpr_meson_compile()
        data = dict()
        data["cexpr_optimized"] = cexpr_optimized
        data["cexpr_original"] = cexpr_original
        data["code_py"] = code_py
        q.save_pickle_obj(data, fn_pickle)
        return cexpr_optimized
    if q.get_id_node() == 0 and not q.does_file_exist(fn_pickle):
        calc_compile_cexpr()
    q.sync_node()
    while not q.does_file_exist(fn_pickle):
        q.displayln(3, f"cache_compiled_cexpr: Node {q.get_id_node()}: waiting for '{fn_pickle}'.")
        time.sleep(0.1)
    cexpr = q.load_pickle_obj(fn_pickle)["cexpr_optimized"]
    module = importlib.import_module((path + "/build/code").replace("/", "."))
    cexpr.function = {
            # cexpr_function(positions_dict, get_prop, is_ama_and_sloppy = False) => val as 1-D np.array
            "cexpr_function" : module.cexpr_function,
            # cexpr_function_get_prop(positions_dict, get_prop) => props
            "cexpr_function_get_prop" : module.cexpr_function_get_prop,
            # cexpr_function_eval(positions_dict, props) => ama_val as AmaVal of 1-D np.array
            "cexpr_function_eval" : module.cexpr_function_eval,
            }
    cexpr.total_sloppy_flops = module.total_sloppy_flops
    return cexpr

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

def benchmark_show_check(check):
    return " ".join([ f"{v:.10E}" for v in check ])

@q.timer
def benchmark_eval_cexpr(cexpr : CExpr, *,
        benchmark_size = 10,
        benchmark_num = 10,
        benchmark_num_ama = 2,
        benchmark_rng_state = None):
    if benchmark_rng_state is None:
        benchmark_rng_state = q.RngState("benchmark_eval_cexpr")
    expr_names = get_cexpr_names(cexpr)
    n_expr = len(expr_names)
    n_pos = len(cexpr.positions)
    prop_dict = {}
    size = q.Coordinate([ 8, 8, 8, 16, ])
    positions = [
            ("point", tuple(benchmark_rng_state.split(f"positions {pos_idx}").c_rand_gen(size).list()),)
            for pos_idx in range(n_pos)
            ]
    set_flavors = set()
    for name, value in cexpr.variables_prop:
        assert name.startswith("V_S_")
        x = value
        assert isinstance(x, Op)
        assert x.otype == "S"
        flavor = x.f
        set_flavors.add(flavor)
    for flavor in set_flavors:
        for pos_src_idx in range(n_pos):
            pos_src = positions[pos_src_idx]
            for pos_snk_idx in range(n_pos):
                pos_snk = positions[pos_snk_idx]
                prop = make_rand_spin_color_matrix(benchmark_rng_state.split(f"prop {flavor} {pos_snk_idx} {pos_src_idx}"))
                prop_ama = make_rand_spin_color_matrix(benchmark_rng_state.split(f"prop ama {flavor} {pos_snk_idx} {pos_src_idx}"))
                prop_dict[(flavor, pos_snk, pos_src,)] = mk_ama_val(prop, pos_src, [ prop, prop_ama, ], [ 0, 1, ], [ 1.0, 0.5, ])
    def mk_pos_dict(k):
        positions_dict = {}
        positions_dict["size"] = size.list()
        idx_list = q.random_permute(list(range(n_pos)), benchmark_rng_state.split(f"pos_dict {k}"))
        for pos, idx in zip(cexpr.positions, idx_list):
            positions_dict[pos] = positions[idx]
        return positions_dict
    positions_dict_list = [ mk_pos_dict(k) for k in range(benchmark_size) ]
    #
    @q.timer
    def get_prop(flavor, pos_snk, pos_src):
        return ama_extract(prop_dict[(flavor, pos_snk, pos_src)], is_sloppy = True)
    @q.timer_verbose
    def benchmark_eval_cexpr_run():
        res_list = []
        for k in range(benchmark_size):
            res = eval_cexpr(cexpr, positions_dict = positions_dict_list[k], get_prop = get_prop)
            res_list.append(res)
        res = np.array(res_list)
        assert res.shape == (benchmark_size, n_expr,)
        return res
    @q.timer
    def get_prop_ama(flavor, pos_snk, pos_src):
        return prop_dict[(flavor, pos_snk, pos_src)]
    @q.timer_verbose
    def benchmark_eval_cexpr_run_with_ama():
        res_list = []
        for k in range(benchmark_size):
            res1 = eval_cexpr(cexpr, positions_dict = positions_dict_list[k], get_prop = get_prop_ama)
            res2 = eval_cexpr(cexpr, positions_dict = positions_dict_list[k], get_prop = get_prop)
            res_ama, res_sloppy = eval_cexpr(cexpr, positions_dict = positions_dict_list[k], get_prop = get_prop_ama, is_ama_and_sloppy = True)
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

def sqr_component(x):
    return x.real * x.real + 1j * x.imag * x.imag

def sqrt_component(x):
    return math.sqrt(x.real) + 1j * math.sqrt(x.imag)

def sqr_component_array(arr):
    return np.array([ sqr_component(x) for x in arr ])

def sqrt_component_array(arr):
    return np.array([ sqrt_component(x) for x in arr ])

def get_mpi_chunk(total_list, *, rng_state = None):
    """
    rng_state has to be the same on all the nodes
    e.g. rng_state = q.RngState("get_mpi_chunk")
    """
    total = len(total_list)
    id_worker = q.get_id_node()
    num_worker = q.get_num_node()
    size_max = (total - 1) // num_worker + 1;
    start = min(id_worker * size_max, total);
    stop = min(start + size_max, total);
    # size = stop - start;
    if rng_state is not None:
        assert isinstance(rng_state, q.RngState)
        total_list = q.random_permute(total_list, rng_state)
    return total_list[start:stop]
