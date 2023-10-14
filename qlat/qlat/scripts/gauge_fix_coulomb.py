# Author: Christoph Lehner 2021
# Author: Luchang Jin 2023

import qlat as q
import sys

if len(sys.argv) == 1:
    q.displayln_info("Usage: gauge-fix-coulomb [ --source source_config ] [ --output gt.field ] [ --guess gt-guess.field ]")
    q.displayln_info("If the output file already exist, it will only check output file.")
    q.displayln_info("Default output file is 'gt.field'.")
    q.displayln_info("Will use CPS gauge fix format is output file end with '.gfix'.")
    sys.exit()

import gpt as g
import qlat_gpt as qg

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
p_guess = g.default.get("--guess", None)
p_output = g.default.get("--output", "gt.field")
p_rng_seed = g.default.get("--random", None)

def load_gf():
    if p_source is None:
        q.displayln_info("Need to provide source file with '--source filename'. Use a sample gauge field for now.")
        total_site = [ 4, 4, 4, 8, ]
        geo = q.Geometry(total_site, 1)
        gf = q.GaugeField(geo)
        rs = q.RngState("seed")
        gf.set_rand(rs.split("gf-init"), 0.5, 10)
    else:
        q.displayln_info(f"Loading gauge field '{p_source}'.")
        if not q.does_file_exist_qar_sync_node(p_source):
            q.displayln_info(f"ERROR: Gauge field file '{p_source}' does not exist.")
            sys.exit(1)
        gf = q.GaugeField()
        gf.load(p_source)
    return gf

def load_gt_guess():
    if p_guess is None:
        q.displayln_info("No guess is provided. Gauge fix from scratch.")
        return None
    elif not q.does_file_exist_qar_sync_node(p_guess):
        q.displayln_info(f"WARNING: Guess file '{p_guess}' does not exist. Gauge fix from scratch.")
        return None
    q.displayln_info(f"Find '{p_guess}'.")
    gt = q.GaugeTransform()
    if p_guess.endswith(".gfix"):
        q.displayln_info(f"Assuming CPS gauge fix format.")
        gt.load_cps(p_guess)
    else:
        q.displayln_info(f"Assuming Qlat field format (double precision, big endian).")
        gt.load_double(p_guess)
    return gt

def load_gt_output():
    assert isinstance(p_output, str)
    if not q.does_file_exist_qar_sync_node(p_output):
        q.displayln_info(f"Output file '{p_output}' does not exist. Will do gauge fix.")
        return None
    q.displayln_info(f"Find '{p_output}'. Will only do checking.")
    gt = q.GaugeTransform()
    if p_output.endswith(".gfix"):
        q.displayln_info(f"Assuming CPS gauge fix format.")
        gt.load_cps(p_output)
    else:
        q.displayln_info(f"Assuming Qlat field format (double precision, big endian).")
        gt.load_double(p_output)
    return gt

def save_gt_output(gt):
    assert isinstance(p_output, str)
    assert not q.does_file_exist_sync_node(p_output)
    q.displayln_info(f"Start to write to '{p_output}'.")
    if p_output.endswith(".gfix"):
        q.displayln_info(f"Assuming CPS gauge fix format.")
        gt.save_cps(p_output)
    else:
        q.displayln_info(f"Assuming Qlat field format (double precision, big endian).")
        gt.save_double(p_output)
    q.sync_node()

gf = load_gf()

gt = load_gt_output()

if gt is None:
    # Do gauge fixing
    gt = load_gt_guess()
    gt = qg.gauge_fix_coulomb(
            gf,
            gt=gt,
            mpi_split=p_mpi_split,
            maxiter_gd=p_maxiter_gd,
            maxiter_cg=p_maxiter_cg,
            maxcycle_cg=p_maxcycle_cg,
            log_every=p_log_every,
            eps=p_eps,
            step=p_step,
            step_gd=p_step_gd,
            rng_seed=p_rng_seed,
            )
    save_gt_output(gt)
else:
    # Only do checking
    qg.check_gauge_fix_coulomb(gf, gt, eps=p_eps)

q.timer_display()

qg.end_with_gpt()

exit()
