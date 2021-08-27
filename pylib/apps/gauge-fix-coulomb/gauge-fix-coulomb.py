#!/usr/bin/env python3
#
# Author: Christoph Lehner 2021
# Author: Luchang Jin 2021

import gpt as g
import qlat as q
import qlat_gpt as qg
import sys

qg.begin_with_gpt()

# Parameters
p_mpi_split = g.default.get_ivec("--mpi_split", [ 1, 1, 1, ], 3)
p_maxiter_cg = g.default.get_int("--maxiter_cg", 200)
p_maxiter_gd = g.default.get_int("--maxiter_gd", 10)
p_maxcycle_cg = g.default.get_int("--maxcycle_cg", 50)
p_log_every = g.default.get_int("--log_every", 1)
p_eps = g.default.get_float("--eps", 1e-12)
p_step = g.default.get_float("--step", 0.3)
p_step_gd = g.default.get_float("--step_gd", 0.1)
p_source = g.default.get("--source", None)
p_output = g.default.get("--output", "gt.field")
p_rng_seed = g.default.get("--random", None)


if p_source is None:
    g.message("Need to provide source file")
    sys.exit(1)

gf = qg.load_gauge_field(p_source)

gt = qg.gauge_fix_coulomb(
        gf,
        mpi_split = p_mpi_split,
        maxiter_gd = p_maxiter_gd,
        maxiter_cg = p_maxiter_cg,
        maxcycle_cg = p_maxcycle_cg,
        log_every = p_log_every,
        eps = p_eps,
        step = p_step,
        step_gd = p_step_gd,
        rng_seed = p_rng_seed,
        )

gt.save_double(p_output)

q.timer_display()

qg.end_with_gpt()
