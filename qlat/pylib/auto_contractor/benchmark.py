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

from auto_contractor.eval import *

import qlat as q
import cqlat

@q.timer
def benchmark_function_1(f, arg, benchmark_size = 1000, benchmark_num = 10, total_flops = 0):
    @q.timer_verbose
    def benchmark_run_10():
        q.acc_timer_flops("py:benchmark_run_10", total_flops * 10 * benchmark_size)
        for k in range(benchmark_size):
            f(arg)
            f(arg)
            f(arg)
            f(arg)
            f(arg)
            f(arg)
            f(arg)
            f(arg)
            f(arg)
            f(arg)
    q.timer_fork(0)
    for i in range(benchmark_num):
        benchmark_run_10()
    q.timer_display()
    q.timer_merge()

@q.timer
def benchmark_function_2(f, arg1, arg2, benchmark_size = 1000, benchmark_num = 10, total_flops = 0):
    @q.timer_verbose
    def benchmark_run_10():
        q.acc_timer_flops("py:benchmark_run_10", total_flops * 10 * benchmark_size)
        for k in range(benchmark_size):
            f(arg1, arg2)
            f(arg1, arg2)
            f(arg1, arg2)
            f(arg1, arg2)
            f(arg1, arg2)
            f(arg1, arg2)
            f(arg1, arg2)
            f(arg1, arg2)
            f(arg1, arg2)
            f(arg1, arg2)
    q.timer_fork(0)
    for i in range(benchmark_num):
        benchmark_run_10()
    q.timer_display()
    q.timer_merge()

@q.timer
def benchmark_function_3(f, arg1, arg2, arg3, benchmark_size = 1000, benchmark_num = 10, total_flops = 0):
    @q.timer_verbose
    def benchmark_run_10():
        q.acc_timer_flops("py:benchmark_run_10", total_flops * 10 * benchmark_size)
        for k in range(benchmark_size):
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
            f(arg1, arg2, arg3)
    q.timer_fork(0)
    for i in range(benchmark_num):
        benchmark_run_10()
    q.timer_display()
    q.timer_merge()

if __name__ == "__main__":
    rng_state = q.RngState("seed")
    mat_sc_1 = make_rand_spin_color_matrix(rng_state)
    mat_sc_2 = make_rand_spin_color_matrix(rng_state)
    mat_sc_3 = make_rand_spin_color_matrix(rng_state)
    mat_s_1 = make_rand_spin_matrix(rng_state)
    mat_s_2 = make_rand_spin_matrix(rng_state)
    mat_s_3 = get_gamma_matrix(1)
    print("sc * sc")
    benchmark_function_2(mat_mul_sc_sc, get_mat(mat_sc_1), get_mat(mat_sc_2), total_flops = 13536)
    print("sc * s")
    benchmark_function_2(mat_mul_sc_s, get_mat(mat_sc_1), get_mat(mat_s_1), total_flops = 4320)
    print("s * sc")
    benchmark_function_2(mat_mul_s_sc, get_mat(mat_s_1), get_mat(mat_sc_1), total_flops = 4320)
    print("s * s")
    benchmark_function_2(mat_mul_s_s, get_mat(mat_s_1), get_mat(mat_s_2), total_flops = 4320)
    print("tr(sc, sc)")
    benchmark_function_2(mat_sc_sc_trace, get_mat(mat_sc_1), get_mat(mat_sc_2), total_flops = 480)
    print("tr(sc)")
    benchmark_function_1(mat_sc_trace, get_mat(mat_sc_1), total_flops = 22)
