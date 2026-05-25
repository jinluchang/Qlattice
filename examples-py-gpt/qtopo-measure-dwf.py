#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g
import numpy as np
import sys

from auto_contractor import (
    mk_scalar,
    mk_scalar5,
    contract_simplify_compile,
    cache_compiled_cexpr,
    benchmark_eval_cexpr,
    get_expr_names,
    eval_cexpr,
)

from qlat_scripts.v1 import (
    set_param,
    get_param,
    is_test,
    run_params,
    run_gf,
    get_load_path,
    get_save_path,
)

usage = f"""
Topological charge measurement based on DWF vector current in the fifth dimension with Qlattice and GPT/Grid
by Luchang Jin
2026/03/21
{""}
{__file__} --usage
# Show this message and exit.
{__file__} --test
# Generate some test data and then perform the measurement.
{""}
{__file__} \
    [--sparse_ratio 32] \
    [--num_of_rand_vol_u1 2] \
    [--ls 16] \
    [--ls_sloppy 12] \
    [--b_plus_c 3.0] \
    [--b_plus_c_sloppy 3.0] \
    [--maxiter_sloppy 50] \
    [--maxiter_exact 100] \
    [--ama_prob 0.1] \
    --gf PATH_GAUGE_FIELD \
    --out PATH_OUTPUT
# Random number seed will be based on PATH_GAUGE_FIELD.
# Multiple input and output paths can be specified (they will be processed the same way).
# The program does not overwrite output files. If the output file already exist, the program will simply display that output file content and skip that input and output file pair and continue.
"""

is_cython = not is_test()
# is_cython = True

rs_sig = q.RngState("rs_sig")

# ------------------------------

@q.timer(is_verbose=True)
def measure_topo_dwf(
    gf,
    info_path,
    *,
    params,
):
    """
    Main measurement routine for topological charge via DWF fifth-dimension vector current.\n
    For each random volume U(1) source, the routine:
      1. Constructs a sparse grid decomposition and solves the DWF propagator on each
         sparse subset using sloppy (and optionally exact AMA-corrected) inversions.
      2. Saves/loads sparse solutions from disk.
      3. Assembles the full propagator and evaluates the tadpole-loop contraction
         expression (qbar gamma5 q and qbar q) at every lattice site.
      4. Sums over time slices and accumulates results across random U(1) sources.\n
    Parameters
    ----------
    gf : q.GaugeField
        The input gauge field configuration.
    info_path : str
        Directory path for saving all output (pickle, JSON, field data, scratch).
    params : dict
        Dictionary of measurement parameters with keys:
          - sparse_ratio (int): coarsening factor for the sparse-grid decomposition.
          - num_of_rand_vol_u1 (int): number of random volume U(1) sources.
          - ls (int): exact Mobius Ls (fifth-dimension length).
          - ls_sloppy (int): sloppy Mobius Ls.
          - b_plus_c (float): exact Mobius b + c parameter.
          - b_plus_c_sloppy (float): sloppy Mobius b + c parameter.
          - maxiter_sloppy (int): max CG iterations for the sloppy solve.
          - maxiter_exact (int): max CG iterations for the exact solve.
          - ama_prob (float): probability of applying the AMA correction.
          - seed (str or int): random seed.\n
    Returns
    -------
    info : dict
        Dictionary with keys:
          - "f_tadpole_loop_sum_arr" (np.ndarray, complex128):
              shape (num_of_rand_vol_u1, nt, 2) where
              [r, t, 0] = topological charge on time-slice t for source r,
              [r, t, 1] = quark condensate on time-slice t for source r.
          - "f_tadpole_loop_imag_sqr_sum_arr" (np.ndarray, float64):
              shape (num_of_rand_vol_u1, nt, 2), imaginary-part squared sums
              used for error estimation.
    """
    q.check_time_limit()
    assert isinstance(gf, q.GaugeField), f"measure_topo_dwf: gf type={type(gf)}"
    assert isinstance(info_path, str), f"measure_topo_dwf: info_path type={type(info_path)}"
    assert isinstance(params, dict), f"measure_topo_dwf: params type={type(params)}"
    lat_shape = tuple(gf.geo.total_site)
    assert len(lat_shape) == 4 and all(s > 0 for s in lat_shape), \
        f"measure_topo_dwf: lat_shape={lat_shape}"
    #
    sparse_ratio = params["sparse_ratio"]
    num_of_rand_vol_u1 = params["num_of_rand_vol_u1"]
    ls = params["ls"]
    ls_sloppy = params["ls_sloppy"]
    b_plus_c = params["b_plus_c"]
    b_plus_c_sloppy = params["b_plus_c_sloppy"]
    maxiter_sloppy = params["maxiter_sloppy"]
    maxiter_exact = params["maxiter_exact"]
    ama_prob = params["ama_prob"]
    seed = params["seed"]
    #
    q.save_json_obj(params, f"{info_path}/params.json")
    #
    rs = q.RngState(seed)
    #
    rand_vol_u1_multiplicity = 12
    #
    mobius_params = {
        "mass": 1.0,
        "M5": 1.8,
        "b": (b_plus_c + 1) / 2,
        "c": (b_plus_c - 1) / 2,
        "Ls": ls,
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
    }
    mobius_params_sloppy = {
        "mass": 1.0,
        "M5": 1.8,
        "b": (b_plus_c_sloppy + 1) / 2,
        "c": (b_plus_c_sloppy - 1) / 2,
        "Ls": ls_sloppy,
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
    }
    #
    mpi_split = g.default.get_ivec("--mpi_split", None, 4)
    if mpi_split is not None:
        n_grouped = g.default.get_int("--grouped", 12)
    else:
        n_grouped = g.default.get_int("--grouped", 1)
    inv = g.algorithms.inverter
    cg_f = inv.split(
        inv.cg({"eps": 1e-8, "maxiter": maxiter_sloppy}), mpi_split=mpi_split
    )
    cg = inv.split(inv.cg({"eps": 1e-8, "maxiter": maxiter_exact}), mpi_split=mpi_split)
    pc = g.qcd.fermion.preconditioner
    #
    slv_5d = inv.defect_correcting(
        inv.mixed_precision(
            inv.preconditioned(pc.eo2_ne(), cg),
            g.single,
            g.double,
        ),
        eps=1e-8,
        maxiter=100,
    )
    slv_5d_f = inv.mixed_precision(
        inv.preconditioned(pc.eo2_ne(), cg_f),
        g.single,
        g.double,
    )
    #
    gf = gf.copy()
    gf.unitarize()
    #
    geo = gf.geo
    total_site = geo.total_site
    #
    q.json_results_append("gf.plaq()", gf.plaq(), 1e-12)
    q.json_results_append("gf.link_trace()", gf.link_trace(), 1e-12)
    #
    gt = q.GaugeTransform(geo)
    gt.set_rand(rs.split("gt-init"), 0.5, 50)
    gt.unitarize()
    gf = gt * gf
    q.json_results_append("after transform gf.plaq()", gf.plaq(), 1e-12)
    q.json_results_append("after transform gf.link_trace()", gf.link_trace(), 1e-12)
    #
    psel_full = q.PointsSelection(geo)
    #
    q.json_results_append(f"{total_site.volume()=}")
    q.json_results_append(f"{psel_full.n_points=}")
    #
    xg_arr_full = psel_full.xg_arr
    #
    psel_list = []
    for idx in range(sparse_ratio):
        sel_arr = mk_sparse_grid(xg_arr_full, sparse_ratio, idx)
        xg_arr = xg_arr_full[sel_arr]
        psel = q.PointsSelection(total_site, xg_arr, "l")
        psel_list.append(psel)
        q.json_results_append(
            f"{sparse_ratio=} {idx=} {len(psel)=} {q.hash_sha256(psel)=}"
        )
    #
    q.json_results_append(f"{len(psel_list)=}")
    #
    gpt_gf = qg.gpt_from_qlat(gf)
    q.json_results_append(
        "g.qcd.gauge.plaquette(gpt_gf)", g.qcd.gauge.plaquette(gpt_gf), 1e-12
    )
    #
    gf1 = qg.qlat_from_gpt(gpt_gf)
    q.json_results_append("gf1.plaq()", gf1.plaq(), 1e-12)
    #
    qm = g.qcd.fermion.mobius(gpt_gf, mobius_params)
    qm_f = g.qcd.fermion.mobius(gpt_gf, mobius_params_sloppy)
    slv_qm = qm.propagator(slv_5d).grouped(n_grouped)
    slv_qm_f = qm_f.propagator(slv_5d_f).grouped(n_grouped)
    inv_qm = qg.InverterGPT(inverter=slv_qm, qtimer=q.Timer("py:slv_qm", True))
    inv_qm_f = qg.InverterGPT(inverter=slv_qm_f, qtimer=q.Timer("py:slv_qm_f", True))
    #
    rs_rand_u1 = rs.split("qtopo-measure-dwf(rand_u1)")
    rs_ama = rs.split("qtopo-measure-dwf(ama)")
    #
    info_list = []
    #
    for rand_vol_u1_idx in range(num_of_rand_vol_u1):
        q.check_time_limit()
        fu1 = q.mk_rand_vol_u1(
            geo, rand_vol_u1_multiplicity, rs_rand_u1.split(f"{rand_vol_u1_idx}")
        )
        fu1[:, 3:6] = fu1[:, 0:3]
        fu1[:, 6:12] = fu1[:, 0:6]
        #
        q.json_results_append(
            "fu1 sig",
            q.get_data_sig_arr(fu1, rs_sig, 3),
            1e-12,
        )
        #
        def mk_path(idx):
            """
            Build the scratch directory path for a given sparse-solve index.\n
            Parameters
            ----------
            idx : int
                Sparse-solve (grid subset) index.\n
            Returns
            -------
            str
                Full scratch path.
            """
            assert isinstance(idx, int), f"mk_path: idx type={type(idx)}"
            return f"{info_path}/scratch/rand_vol_u1_idx-{rand_vol_u1_idx}/sparse_solve_idx-{idx}"
        #
        def check(idx):
            """
            Check whether scratch data (psel + sparse propagator) already exists on disk.\n
            Parameters
            ----------
            idx : int
                Sparse-solve index.\n
            Returns
            -------
            bool
                True if both ``psel.lati`` and ``sp_prop_sol.lat`` exist.
            """
            assert isinstance(idx, int), f"check: idx type={type(idx)}"
            if not q.does_file_exist_qar_sync_node(f"{mk_path(idx)}/psel.lati"):
                return False
            if not q.does_file_exist_qar_sync_node(f"{mk_path(idx)}/sp_prop_sol.lat"):
                return False
            return True
        #
        @q.timer(is_verbose=True)
        def save(idx, sp_prop_sol):
            """
            Shuffle the sparse propagator to the root node and save it along with its point selection.\n
            Only MPI rank 0 actually writes to disk; other ranks assert that their
            shuffled data is empty.\n
            Parameters
            ----------
            idx : int
                Sparse-solve index (used to build the output path).
            sp_prop_sol : q.PselProp
                Sparse propagator solution in the local (``"l"``) point distribution.
            """
            assert isinstance(idx, int), f"save: idx type={type(idx)}"
            assert isinstance(sp_prop_sol, q.PselProp), f"save: sp_prop_sol type={type(sp_prop_sol)}"
            root = 0
            psel = psel_list[idx]
            assert sp_prop_sol.n_points == psel.n_points
            ssp = q.SelectedShufflePlan("g_from_l", psel, root, geo)
            psel_g = ssp.psel_recv_list[0]
            sp_prop_sol_g = ssp.shuffle_sp(q.PselProp, sp_prop_sol)
            if root == q.get_id_node():
                assert sp_prop_sol_g.n_points == psel_g.n_points
                psel_g.save(
                    f"{mk_path(idx)}/psel.lati",
                    is_sync_node=False,
                )
                sp_prop_sol_g.save(
                    f"{mk_path(idx)}/sp_prop_sol.lat",
                    is_sync_node=False,
                )
            else:
                assert len(psel_g) == 0
                assert sp_prop_sol_g.n_points == 0
        #
        @q.timer(is_verbose=True)
        def load(idx_list):
            """
            Load sparse propagator solutions from disk and redistribute across MPI ranks.\n
            Each index is assigned a root rank via round-robin over the shuffle node list.
            After reading on the respective root, the data is reverse-shuffled back to
            the local (``"l"``) point distribution.\n
            Parameters
            ----------
            idx_list : list of int
                Sparse-solve indices to load.\n
            Returns
            -------
            sp_prop_sol_il : list of q.PselProp
                List of loaded sparse propagator solutions, one per index.
            """
            assert isinstance(idx_list, list), f"load: idx_list type={type(idx_list)}"
            for i_idx in idx_list:
                assert isinstance(i_idx, int), f"load: idx_list element type={type(i_idx)}"
            id_node_list_for_shuffle = q.get_id_node_list_for_shuffle()
            root_il = [
                id_node_list_for_shuffle[i % len(id_node_list_for_shuffle)]
                for i, idx in enumerate(idx_list)
            ]
            geo_il = [geo for idx in idx_list]
            psel_il = [psel_list[idx] for idx in idx_list]
            ssp = q.SelectedShufflePlan("g_from_l", psel_il, root_il, geo_il)
            psel_g_il = ssp.psel_recv_list
            sp_prop_sol_g_il = [
                q.PselProp(psel_g_il[i]) for i, idx in enumerate(idx_list)
            ]
            for i, idx in enumerate(idx_list):
                if root_il[i] == q.get_id_node():
                    psel_g_l = q.PointsSelection()
                    psel_g_l.load(
                        f"{mk_path(idx)}/psel.lati",
                        is_sync_node=False,
                    )
                    assert psel_g_il[i] == psel_g_l
                    sp_prop_sol_g_il[i].load(
                        f"{mk_path(idx)}/sp_prop_sol.lat",
                        is_sync_node=False,
                    )
            sp_prop_sol_il = ssp.shuffle_sp_list(
                q.PselProp,
                sp_prop_sol_g_il,
                is_reverse=True,
            )
            assert len(sp_prop_sol_il) == len(idx_list)
            assert isinstance(sp_prop_sol_il, list), f"load: sp_prop_sol_il type={type(sp_prop_sol_il)}"
            for sp in sp_prop_sol_il:
                assert isinstance(sp, q.PselProp), f"load: sp_prop_sol_il element type={type(sp)}"
            return sp_prop_sol_il
        #
        @q.timer
        def benchmark(sp_prop_sol):
            """
            Log quality metrics (signal/error ratio) for a sparse propagator solution.\n
            Uses gamma_5 Hermiticity to separate signal (gamma_5-hermitian part) from
            noise (anti-hermitian part). Reports:
              - prop_size: Frobenius-norm signal, error, and error/signal ratio.
              - topo_sum: trace(gamma_5 * prop), scaled by sqrt(volume).
              - quark_condensate_avg: trace(prop) averaged over points.\n
            Parameters
            ----------
            sp_prop_sol : q.PselProp
                Sparse propagator solution to benchmark.
            """
            assert isinstance(sp_prop_sol, q.PselProp), f"benchmark: sp_prop_sol type={type(sp_prop_sol)}"
            fname = q.get_fname()
            n_points = q.glb_sum(sp_prop_sol.n_points)
            psel = sp_prop_sol.psel
            total_site = psel.total_site
            volume = total_site.volume()
            gamma_5_s = q.get_gamma_matrix(5)[:]
            gamma_5_sc = np.zeros(
                (
                    4,
                    3,
                    4,
                    3,
                ),
                dtype=np.complex128,
            )
            gamma_5_sc[:, np.arange(3), :, np.arange(3)] = gamma_5_s[None, :, :]
            gamma_5 = gamma_5_sc.reshape(12, 12).copy()
            #
            arr = sp_prop_sol[:, 0, :, :]
            arr_t = gamma_5 @ arr.conj().swapaxes(-1, -2) @ gamma_5
            arr_avg = 0.5 * (arr + arr_t)
            arr_diff = 0.5 * (arr - arr_t)
            val = np.sqrt(q.glb_sum(q.qnorm(arr_avg)) / n_points)
            err = np.sqrt(q.glb_sum(q.qnorm(arr_diff)) / n_points)
            q.json_results_append(
                f"{fname}: prop_size val,err,err/val",
                np.array(
                    [
                        val,
                        err,
                        err / val,
                    ]
                ),
                1e-5,
            )
            topo_arr = np.trace(gamma_5 @ arr, axis1=-1, axis2=-2)
            arr_avg = topo_arr.real
            arr_diff = topo_arr.imag
            val = np.sqrt(q.glb_sum(q.qnorm(arr_avg)) / n_points * volume)
            err = np.sqrt(q.glb_sum(q.qnorm(arr_diff)) / n_points * volume)
            q.json_results_append(
                f"{fname}: topo_sum val,err,err/val",
                np.array(
                    [
                        val,
                        err,
                        err / val,
                    ]
                ),
                1e-5,
            )
            quark_condensate_arr = np.trace(arr, axis1=-1, axis2=-2)
            arr_avg = quark_condensate_arr.real
            arr_diff = quark_condensate_arr.imag
            val = np.sqrt(q.glb_sum(q.qnorm(arr_avg)) / n_points)
            err = np.sqrt(q.glb_sum(q.qnorm(arr_diff)) / n_points / volume)
            q.json_results_append(
                f"{fname}: quark_condensate_avg val,err,err/val",
                np.array(
                    [
                        val,
                        err,
                        err / val,
                    ]
                ),
                1e-5,
            )
        #
        for idx, psel in enumerate(psel_list):
            q.check_time_limit()
            q.json_results_append(f"sparse_solve: {idx + 1}/{len(psel_list)}")
            if check(idx):
                continue
            sp_prop_sol = sparse_solve(idx, psel, fu1, inv_qm_f)
            benchmark(sp_prop_sol)
            rs_ama_idx = rs_ama.split(f"{rand_vol_u1_idx} {idx}")
            r = rs_ama_idx.u_rand_gen()
            if r <= ama_prob:
                sp_prop_sol_ama = sparse_solve(idx, psel, fu1, inv_qm)
                sp_prop_sol_ama -= sp_prop_sol
                q.json_results_append(
                    f"sp_prop_sol_ama ratio of sqrt(qnorm) ({idx + 1}/{len(psel_list)}) ({ama_prob=:.4f})",
                    np.sqrt(
                        q.glb_sum(q.qnorm(sp_prop_sol_ama))
                        / q.glb_sum(q.qnorm(sp_prop_sol))
                    ).item(),
                    1e-5,
                )
                sp_prop_sol_ama *= 1 / ama_prob
                sp_prop_sol += sp_prop_sol_ama
                benchmark(sp_prop_sol)
            save(idx, sp_prop_sol)
            sp_prop_sol_il = load(
                [
                    idx,
                ]
            )
            assert np.all(
                q.get_data_sig_arr(sp_prop_sol, rs_sig, 3)
                == q.get_data_sig_arr(sp_prop_sol_il[0], rs_sig, 3)
            )
        #
        sp_prop_sol_list = load([idx for idx in range(len(psel_list))])
        #
        prop_sol = q.Prop(geo)
        prop_sol.set_zero()
        for sp_prop_sol in sp_prop_sol_list:
            prop_sol @= sp_prop_sol
        q.json_results_append(
            "prop_sol sig",
            q.get_data_sig_arr(prop_sol, rs_sig, 3),
            1e-7,
        )
        prop_sol.save_float_from_double(
            f"{info_path}/field/rand_vol_u1_idx-{rand_vol_u1_idx}/prop_sol.field",
        )
        #
        def get_prop(flavor, p_snk, p_src):
            """
            Callback for the contraction expression evaluator: return the propagator element.\n
            Only the ``"c"`` (connected) flavor is supported.  Sink and source must both
            be ``"point-snk"`` at the same coordinate.\n
            Parameters
            ----------
            flavor : str
                Flavor channel (must be ``"c"``).
            p_snk : tuple(str, tuple)
                Sink specification: ``("point-snk", (x, y, z, t))``.
            p_src : tuple(str, tuple)
                Source specification: ``("point-snk", (x, y, z, t))``.\n
            Returns
            -------
            WilsonMatrix
                12x12 Dirac-colour propagator matrix at the given site.
            """
            assert isinstance(flavor, str), f"get_prop: flavor type={type(flavor)}"
            assert flavor == "c", f"get_prop: flavor={flavor} != 'c'"
            assert isinstance(p_snk, tuple) and isinstance(p_src, tuple), \
                f"get_prop: p_snk type={type(p_snk)} p_src type={type(p_src)}"
            type_snk, pos_snk = p_snk
            type_src, pos_src = p_src
            assert type_snk == "point-snk", f"get_prop: type_snk={type_snk}"
            assert type_src == "point-snk", f"get_prop: type_src={type_src}"
            pos_snk = q.Coordinate(pos_snk)
            pos_src = q.Coordinate(pos_src)
            assert pos_snk == pos_src, f"get_prop: pos_snk={pos_snk} != pos_src={pos_src}"
            index = geo.index_from_g_coordinate(pos_snk)
            val = prop_sol.get_elem_wm(index)
            assert isinstance(val, q.WilsonMatrix), \
                f"get_prop: type(val)={type(val)}"
            return val
        #
        cexpr = get_cexpr_tadpole_loop()
        expr_names = get_expr_names(cexpr)
        chunk_list = q.get_chunk_list(
            xg_arr_full,
            chunk_number=max(1, q.get_q_num_mp_processes()),
        )
        #
        @q.timer(is_verbose=True)
        def eval(chunk_idx):
            """
            Evaluate the tadpole-loop contraction over a chunk of grid coordinates.\n
            For each coordinate in the chunk, builds a position dictionary and calls
            ``eval_cexpr`` with the tadpole-loop compiled expression and the
            ``get_prop`` callback.\n
            Parameters
            ----------
            chunk_idx : int
                Index into ``chunk_list``, defining the subset of coordinates to process.\n
            Returns
            -------
            val_arr : np.ndarray, complex128
                Shape ``(len(chunk), len(expr_names))``.
                Column 0 = ``qbar gamma5 q``, column 1 = ``qbar q`` at each site.
            """
            assert isinstance(chunk_idx, int), f"eval: chunk_idx type={type(chunk_idx)}"
            chunk = chunk_list[chunk_idx]
            val_arr = np.zeros(
                (
                    len(chunk),
                    len(expr_names),
                ),
                dtype=np.complex128,
            )
            for idx, xg in enumerate(chunk):
                pd = {"x_1": ("point-snk", tuple(xg))}
                val_arr[idx] = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            expected_shape = (len(chunk), len(expr_names))
            assert val_arr.shape == expected_shape, f"eval: val_arr.shape={val_arr.shape} != {expected_shape}"
            assert val_arr.dtype == np.complex128
            return val_arr
        #
        val_arr_list = q.parallel_map(eval, range(len(chunk_list)))
        #
        f_tadpole_loop = q.FieldComplexD(geo, 2)
        f_tadpole_loop.set_zero()
        for idx, chunk in enumerate(chunk_list):
            f_tadpole_loop.set_elems_xg(chunk, val_arr_list[idx])
        #
        q.json_results_append(
            "f_tadpole_loop sig",
            q.get_data_sig_arr(f_tadpole_loop, rs_sig, 3),
            1e-7,
        )
        f_tadpole_loop.save_double(
            f"{info_path}/field/rand_vol_u1_idx-{rand_vol_u1_idx}/f_tadpole_loop.field",
        )
        #
        f_tadpole_loop_sum = f_tadpole_loop.glb_sum_tslice()[:].copy()
        q.json_results_append(
            "f_tadpole_loop_sum sig",
            q.get_data_sig_arr(f_tadpole_loop_sum, rs_sig, 3),
            1e-7,
        )
        f_tadpole_loop_imag_sqr = q.FieldRealD(geo, 2)
        f_tadpole_loop_imag_sqr[:] = f_tadpole_loop[:].imag ** 2
        f_tadpole_loop_imag_sqr_sum = f_tadpole_loop_imag_sqr.glb_sum_tslice()[:].copy()
        #
        info = dict()
        info["f_tadpole_loop_sum"] = f_tadpole_loop_sum
        info["f_tadpole_loop_imag_sqr_sum"] = f_tadpole_loop_imag_sqr_sum
        #
        q.save_pickle_obj(
            info,
            f"{info_path}/pickle/rand_vol_u1_idx-{rand_vol_u1_idx}/info.pickle",
        )
        #
        q.save_json_obj(
            dict(
                total=f_tadpole_loop_sum[:, 0].real.sum().item(),
                total_err=np.sqrt(f_tadpole_loop_imag_sqr_sum[:, 0].sum()).item(),
                tslice_sum=f_tadpole_loop_sum[:, 0].real,
                tslice_sum_err=np.sqrt(f_tadpole_loop_imag_sqr_sum[:, 0].real),
            ),
            f"{info_path}/info/rand_vol_u1_idx-{rand_vol_u1_idx}/topo.json",
        )
        q.save_json_obj(
            dict(
                avg=f_tadpole_loop_sum[:, 1].real.sum().item() / geo.total_volume,
                avg_err=np.sqrt(f_tadpole_loop_imag_sqr_sum[:, 1].sum()).item()
                / geo.total_volume,
                tslice_avg=f_tadpole_loop_sum[:, 1].real / geo.spatial_volume,
                tslice_avg_err=np.sqrt(f_tadpole_loop_imag_sqr_sum[:, 1].real)
                / geo.spatial_volume,
            ),
            f"{info_path}/info/rand_vol_u1_idx-{rand_vol_u1_idx}/quark_condensate.json",
        )
        #
        q.json_results_append(
            "f_tadpole_loop_sum topo sum and err",
            np.array(
                [
                    f_tadpole_loop_sum[:, 0].sum(-1),
                    np.sqrt(f_tadpole_loop_imag_sqr_sum[:, 0].sum()).item(),
                ]
            ),
            1e-7,
        )
        q.json_results_append(
            "f_tadpole_loop_sum quark condensate avg and err",
            np.array(
                [
                    f_tadpole_loop_sum[:, 1].sum(-1) / geo.total_volume,
                    np.sqrt(f_tadpole_loop_imag_sqr_sum[:, 1].sum()).item()
                    / geo.total_volume,
                ]
            ),
            1e-7,
        )
        #
        info_list.append(info)
    #
    f_tadpole_loop_sum_arr = np.array(
        [info["f_tadpole_loop_sum"] for info in info_list],
        dtype=np.complex128,
    )
    f_tadpole_loop_imag_sqr_sum_arr = np.array(
        [info["f_tadpole_loop_imag_sqr_sum"] for info in info_list],
        dtype=np.float64,
    )
    #
    info = dict()
    info["f_tadpole_loop_sum_arr"] = f_tadpole_loop_sum_arr
    info["f_tadpole_loop_imag_sqr_sum_arr"] = f_tadpole_loop_imag_sqr_sum_arr
    #
    q.save_pickle_obj(
        info,
        f"{info_path}/pickle/info.pickle",
    )
    #
    q.save_json_obj(
        dict(
            total=f_tadpole_loop_sum_arr[:, :, 0].real.sum(-1),
            total_err=np.sqrt(f_tadpole_loop_imag_sqr_sum_arr[:, :, 0].sum(-1)),
            tslice_sum=f_tadpole_loop_sum_arr[:, :, 0].real,
            tslice_sum_err=np.sqrt(f_tadpole_loop_imag_sqr_sum_arr[:, :, 0]),
        ),
        f"{info_path}/info/topo.json",
    )
    q.save_json_obj(
        dict(
            avg=f_tadpole_loop_sum_arr[:, :, 1].real.sum(-1) / geo.total_volume,
            avg_err=np.sqrt(f_tadpole_loop_imag_sqr_sum_arr[:, :, 1].sum(-1))
            / geo.total_volume,
            tslice_avg=f_tadpole_loop_sum_arr[:, :, 1].real / geo.spatial_volume,
            tslice_avg_err=np.sqrt(f_tadpole_loop_imag_sqr_sum_arr[:, :, 1])
            / geo.spatial_volume,
        ),
        f"{info_path}/info/quark_condensate.json",
    )
    #
    q.json_results_append(
        "f_tadpole_loop_sum_arr sig",
        q.get_data_sig_arr(f_tadpole_loop_sum_arr, rs_sig, 3),
        1e-7,
    )
    q.json_results_append(
        "f_tadpole_loop_sum_arr topo sum and err",
        np.array(
            [
                f_tadpole_loop_sum_arr[:, :, 0].sum(-1),
                np.sqrt(f_tadpole_loop_imag_sqr_sum_arr[:, :, 0].sum(-1)),
            ]
        ),
        1e-7,
    )
    q.json_results_append(
        "f_tadpole_loop_sum_arr quark condensate avg and err",
        np.array(
            [
                f_tadpole_loop_sum_arr[:, :, 1].sum(-1) / geo.total_volume,
                np.sqrt(f_tadpole_loop_imag_sqr_sum_arr[:, :, 1].sum(-1))
                / geo.total_volume,
            ]
        ),
        1e-7,
    )
    #
    q.qtouch_info(f"{info_path}/checkpoint.txt")
    q.qremove_all_info(f"{info_path}/scratch")
    #
    nt = geo.total_site[3]
    assert f_tadpole_loop_sum_arr.shape == (num_of_rand_vol_u1, nt, 2), \
        f"Expected shape ({num_of_rand_vol_u1}, {nt}, 2), got {f_tadpole_loop_sum_arr.shape}"
    assert f_tadpole_loop_imag_sqr_sum_arr.shape == (num_of_rand_vol_u1, nt, 2), \
        f"Expected shape ({num_of_rand_vol_u1}, {nt}, 2), got {f_tadpole_loop_imag_sqr_sum_arr.shape}"
    assert f_tadpole_loop_sum_arr.dtype == np.complex128
    assert f_tadpole_loop_imag_sqr_sum_arr.dtype == np.float64
    #
    assert isinstance(info, dict), f"measure_topo_dwf: info type={type(info)}"
    assert "f_tadpole_loop_sum_arr" in info
    assert isinstance(info["f_tadpole_loop_sum_arr"], np.ndarray)
    assert info["f_tadpole_loop_sum_arr"].dtype == np.complex128
    assert "f_tadpole_loop_imag_sqr_sum_arr" in info
    assert isinstance(info["f_tadpole_loop_imag_sqr_sum_arr"], np.ndarray)
    assert info["f_tadpole_loop_imag_sqr_sum_arr"].dtype == np.float64
    #
    return info

# ------------------------------

@q.timer
def get_cexpr_tadpole_loop():
    """
    Return the compiled contraction expression for the tadpole-loop.\n
    The expression contains two operators evaluated at a single point ``x_1``:
      - ``qbar gamma5 q``  (topological charge density)
      - ``qbar q``         (quark condensate)\n
    The expression is cached to disk under ``cache/auto_contract_cexpr/``
    so that recompilation is avoided on subsequent runs.\n
    Returns
    -------
    cexpr : auto_contractor compiled expression
        Compiled contraction expression object ready for evaluation.
    """
    fn_base = "cache/auto_contract_cexpr/get_cexpr_tadpole_loop"
    #
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[(("x_1", "x_1"), 1)] = "Type1"
        exprs = [
            mk_scalar5("c", "c", "x_1") + "qbar g5 q",
            mk_scalar("c", "c", "x_1") + "qbar q",
        ]
        cexpr = contract_simplify_compile(
            *exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict
        )
        return cexpr
    #
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_verbose=True)
def get_all_cexpr():
    """
    Benchmark the evaluation performance of the tadpole-loop contraction expression.\n
    Results are logged via ``benchmark_eval_cexpr`` for diagnostic and tuning purposes.
    """
    benchmark_eval_cexpr(get_cexpr_tadpole_loop())

@q.timer(is_verbose=True)
def mk_sparse_grid(xg_arr, sparse_ratio, idx):
    """
    Build a boolean selection mask for a sparse sub-grid of the full coordinate array.\n
    The lattice is divided into ``sparse_ratio`` interleaved sub-grids; ``idx`` selects
    which sub-grid to return.  Supported ratios: 1, 2, 4, 8, 16, 32, 64, 81, 128, 162,
    256, 512.\n
    For power-of-two ratios the sub-grids are hyper-rectangular checkerboard patterns;
    for ratio 81 / 162 / 256 / 512 they are based on mod-3 / mod-4 block decomposition.\n
    Parameters
    ----------
    xg_arr : np.ndarray, int
        Shape ``(N, 4)`` array of global grid coordinates ``(x, y, z, t)``.
    sparse_ratio : int
        Total number of sparse sub-grids.
    idx : int
        Index of the desired sub-grid, ``0 <= idx < sparse_ratio``.\n
    Returns
    -------
    sel2_arr : np.ndarray, bool
        Boolean mask of shape ``(N,)``, ``True`` for points belonging to the selected
        sub-grid.
    """
    assert isinstance(xg_arr, np.ndarray), f"mk_sparse_grid: xg_arr type={type(xg_arr)}"
    assert isinstance(sparse_ratio, int), f"mk_sparse_grid: sparse_ratio type={type(sparse_ratio)}"
    assert isinstance(idx, int), f"mk_sparse_grid: idx type={type(idx)}"
    n_points = len(xg_arr)
    expected_xg_shape = (n_points, 4)
    expected_sel_shape = (n_points,)
    assert xg_arr.shape == expected_xg_shape, \
        f"mk_sparse_grid: xg_arr.shape={xg_arr.shape} != {expected_xg_shape}"
    assert xg_arr.dtype in (np.int32, np.int64), f"mk_sparse_grid: xg_arr.dtype={xg_arr.dtype}"
    if sparse_ratio == 1:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        sel_arr &= (xg_arr[:, 0] + xg_arr[:, 1] + xg_arr[:, 2] + xg_arr[:, 3]) % 1 == 0
        sel2_arr = sel_arr.copy()
    elif sparse_ratio == 2:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        sel_arr &= (
            xg_arr[:, 0] + xg_arr[:, 1] + xg_arr[:, 2] + xg_arr[:, 3] + idx0
        ) % 2 == 0
        sel2_arr = sel_arr.copy()
    elif sparse_ratio == 4:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        idx1 = (idx // 2) % 2
        sel_arr &= (xg_arr[:, 0] + xg_arr[:, 1] + idx0) % 2 == 0
        sel_arr &= (xg_arr[:, 2] + xg_arr[:, 3] + idx1) % 2 == 0
        sel2_arr = sel_arr.copy()
        sel2_arr &= (
            xg_arr[:, 0] + xg_arr[:, 1] + xg_arr[:, 2] + xg_arr[:, 3] + idx0 + idx1
        ) % 2 == 0
    elif sparse_ratio == 8:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        idx1 = (idx // 2) % 2
        idx2 = (idx // 4) % 2
        sel_arr &= (xg_arr[:, 0] + xg_arr[:, 1] + idx0) % 2 == 0
        sel_arr &= (xg_arr[:, 2] + xg_arr[:, 3] + idx1) % 2 == 0
        sel_arr &= (xg_arr[:, 0] + xg_arr[:, 2] + idx2) % 2 == 0
        sel2_arr = sel_arr.copy()
        sel2_arr &= (xg_arr[:, 1] - xg_arr[:, 2] + idx0 - idx2) % 2 == 0
        sel2_arr &= (xg_arr[:, 0] - xg_arr[:, 3] + idx2 - idx1) % 2 == 0
        sel2_arr &= (
            xg_arr[:, 0] + xg_arr[:, 1] + xg_arr[:, 2] + xg_arr[:, 3] + idx0 + idx1
        ) % 2 == 0
    elif sparse_ratio == 16:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        idx1 = (idx // 2) % 2
        idx2 = (idx // 4) % 2
        idx3 = (idx // 8) % 2
        sel_arr &= (xg_arr[:, 0] + idx0) % 2 == 0
        sel_arr &= (xg_arr[:, 1] + idx1) % 2 == 0
        sel_arr &= (xg_arr[:, 2] + idx2) % 2 == 0
        sel_arr &= (xg_arr[:, 3] + idx3) % 2 == 0
        sel2_arr = sel_arr.copy()
    elif sparse_ratio == 32:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        idx1 = (idx // 2) % 2
        idx2 = (idx // 4) % 2
        idx3 = (idx // 8) % 2
        idx4 = (idx // 16) % 2
        sel_arr &= (xg_arr[:, 0] + idx0) % 2 == 0
        sel_arr &= (xg_arr[:, 1] + idx1) % 2 == 0
        sel_arr &= (xg_arr[:, 2] + idx2) % 2 == 0
        sel_arr &= (xg_arr[:, 3] + idx3) % 2 == 0
        sel_arr &= (
            (xg_arr[:, 0] + idx0) // 2
            + (xg_arr[:, 1] + idx1) // 2
            + (xg_arr[:, 2] + idx2) // 2
            + (xg_arr[:, 3] + idx3) // 2
            + idx4
        ) % 2 == 0
        sel2_arr = sel_arr.copy()
    elif sparse_ratio == 64:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        idx1 = (idx // 2) % 2
        idx2 = (idx // 4) % 2
        idx3 = (idx // 8) % 2
        idx4 = (idx // 16) % 2
        idx5 = (idx // 32) % 2
        sel_arr &= (xg_arr[:, 0] + idx0) % 2 == 0
        sel_arr &= (xg_arr[:, 1] + idx1) % 2 == 0
        sel_arr &= (xg_arr[:, 2] + idx2) % 2 == 0
        sel_arr &= (xg_arr[:, 3] + idx3) % 2 == 0
        sel_arr &= (
            (xg_arr[:, 0] + idx0) // 2 + (xg_arr[:, 1] + idx1) // 2 + idx4
        ) % 2 == 0
        sel_arr &= (
            (xg_arr[:, 2] + idx2) // 2 + (xg_arr[:, 3] + idx3) // 2 + idx5
        ) % 2 == 0
        sel2_arr = sel_arr.copy()
        sel2_arr &= (
            (xg_arr[:, 0] + idx0) // 2
            + (xg_arr[:, 1] + idx1) // 2
            + (xg_arr[:, 2] + idx2) // 2
            + (xg_arr[:, 3] + idx3) // 2
            + idx4
            + idx5
        ) % 2 == 0
    elif sparse_ratio == 128:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        idx1 = (idx // 2) % 2
        idx2 = (idx // 4) % 2
        idx3 = (idx // 8) % 2
        idx4 = (idx // 16) % 2
        idx5 = (idx // 32) % 2
        idx6 = (idx // 64) % 2
        sel_arr &= (xg_arr[:, 0] + idx0) % 2 == 0
        sel_arr &= (xg_arr[:, 1] + idx1) % 2 == 0
        sel_arr &= (xg_arr[:, 2] + idx2) % 2 == 0
        sel_arr &= (xg_arr[:, 3] + idx3) % 2 == 0
        sel_arr &= (
            (xg_arr[:, 0] + idx0) // 2 + (xg_arr[:, 1] + idx1) // 2 + idx4
        ) % 2 == 0
        sel_arr &= (
            (xg_arr[:, 2] + idx2) // 2 + (xg_arr[:, 3] + idx3) // 2 + idx5
        ) % 2 == 0
        sel_arr &= (
            (xg_arr[:, 0] + idx0) // 2 + (xg_arr[:, 2] + idx2) // 2 + idx6
        ) % 2 == 0
        sel2_arr = sel_arr.copy()
        sel2_arr &= (
            (xg_arr[:, 1] + idx1) // 2 - (xg_arr[:, 2] + idx2) // 2 + idx4 - idx6
        ) % 2 == 0
        sel2_arr &= (
            (xg_arr[:, 0] + idx0) // 2 - (xg_arr[:, 3] + idx3) // 2 + idx6 - idx5
        ) % 2 == 0
        sel2_arr &= (
            (xg_arr[:, 0] + idx0) // 2
            + (xg_arr[:, 1] + idx1) // 2
            + (xg_arr[:, 2] + idx2) // 2
            + (xg_arr[:, 3] + idx3) // 2
            + idx4
            + idx5
        ) % 2 == 0
    elif sparse_ratio == 81:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 3
        idx1 = (idx // 3) % 3
        idx2 = (idx // 9) % 3
        idx3 = (idx // 27) % 3
        sel_arr &= (xg_arr[:, 0] + idx0) % 3 == 0
        sel_arr &= (xg_arr[:, 1] + idx1) % 3 == 0
        sel_arr &= (xg_arr[:, 2] + idx2) % 3 == 0
        sel_arr &= (xg_arr[:, 3] + idx3) % 3 == 0
        sel2_arr = sel_arr.copy()
    elif sparse_ratio == 162:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 3
        idx1 = (idx // 3) % 3
        idx2 = (idx // 9) % 3
        idx3 = (idx // 27) % 3
        idx4 = (idx // 81) % 2
        sel_arr &= (xg_arr[:, 0] + idx0) % 3 == 0
        sel_arr &= (xg_arr[:, 1] + idx1) % 3 == 0
        sel_arr &= (xg_arr[:, 2] + idx2) % 3 == 0
        sel_arr &= (xg_arr[:, 3] + idx3) % 3 == 0
        sel_arr &= (
            (xg_arr[:, 0] + idx0) // 3
            + (xg_arr[:, 1] + idx1) // 3
            + (xg_arr[:, 2] + idx2) // 3
            + (xg_arr[:, 3] + idx3) // 3
            + idx4
        ) % 2 == 0
        sel2_arr = sel_arr.copy()
    elif sparse_ratio == 256:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 4
        idx1 = (idx // 4) % 4
        idx2 = (idx // 16) % 4
        idx3 = (idx // 64) % 4
        sel_arr &= (xg_arr[:, 0] + idx0) % 4 == 0
        sel_arr &= (xg_arr[:, 1] + idx1) % 4 == 0
        sel_arr &= (xg_arr[:, 2] + idx2) % 4 == 0
        sel_arr &= (xg_arr[:, 3] + idx3) % 4 == 0
        sel2_arr = sel_arr.copy()
    elif sparse_ratio == 512:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 4
        idx1 = (idx // 4) % 4
        idx2 = (idx // 16) % 4
        idx3 = (idx // 64) % 4
        idx4 = (idx // 256) % 2
        sel_arr &= (xg_arr[:, 0] + idx0) % 4 == 0
        sel_arr &= (xg_arr[:, 1] + idx1) % 4 == 0
        sel_arr &= (xg_arr[:, 2] + idx2) % 4 == 0
        sel_arr &= (xg_arr[:, 3] + idx3) % 4 == 0
        sel_arr &= (
            (xg_arr[:, 0] + idx0) // 4
            + (xg_arr[:, 1] + idx1) // 4
            + (xg_arr[:, 2] + idx2) // 4
            + (xg_arr[:, 3] + idx3) // 4
            + idx4
        ) % 2 == 0
        sel2_arr = sel_arr.copy()
    else:
        raise Exception(f"{sparse_ratio=}")
    assert isinstance(sel2_arr, np.ndarray), f"mk_sparse_grid: sel2_arr type={type(sel2_arr)}"
    assert sel2_arr.dtype == np.bool_, f"mk_sparse_grid: sel2_arr.dtype={sel2_arr.dtype}"
    assert np.all(sel2_arr == sel_arr)
    assert sel2_arr.shape == expected_sel_shape, \
        f"mk_sparse_grid: sel2_arr.shape={sel2_arr.shape} != {expected_sel_shape}"
    return sel2_arr

@q.timer(is_verbose=True)
def sparse_solve(idx, psel, fu1, inverter):
    """
    Solve the DWF propagator on a sparse subset of the lattice with a random U(1) source.\n
    The source is constructed by projecting the volume U(1) field onto the sparse
    point selection, then solving with the provided inverter.  The solution is
    multiplied by the complex conjugate of the random phase to undo the initial
    phase factor.\n
    Parameters
    ----------
    idx : int
        Sparse-solve index (used for logging only).
    psel : q.PointsSelection
        Sparse point selection defining the source locations.
    fu1 : q.FieldComplexD
        Volume random U(1) field (multiplicity 12).
    inverter : qg.InverterGPT
        GPT inverter for the DWF operator (sloppy or exact).\n
    Returns
    -------
    sp_prop_sol : q.PselProp
        Sparse propagator solution with the conjugate random phase already multiplied.
    """
    assert isinstance(idx, int), f"sparse_solve: idx type={type(idx)}"
    assert isinstance(psel, q.PointsSelection), f"sparse_solve: psel type={type(psel)}"
    assert isinstance(fu1, q.FieldComplexD), f"sparse_solve: fu1 type={type(fu1)}"
    assert isinstance(inverter, qg.InverterGPT), f"sparse_solve: inverter type={type(inverter)}"
    assert fu1.multiplicity == 12, f"sparse_solve: fu1.multiplicity={fu1.multiplicity} != 12"
    geo = fu1.geo
    sp_fu1 = q.SelectedPointsComplexD(psel, fu1.multiplicity)
    sp_fu1 @= fu1
    sp_prop_src = q.PselProp(psel)
    sp_prop_src.set_zero()
    sp_prop_src[:, :, np.arange(12), np.arange(12)] = sp_fu1[:, None, np.arange(12)]
    if is_test():
        q.json_results_append(
            f"sp_prop_src {idx}", q.get_data_sig_arr(sp_prop_src, rs_sig, 3), 1e-12
        )
    sparse_grid_prop_src = q.Prop(geo)
    sparse_grid_prop_src.set_zero()
    sparse_grid_prop_src @= sp_prop_src
    sparse_grid_prop_sol = inverter * sparse_grid_prop_src
    sp_prop_sol = q.PselProp(psel)
    sp_prop_sol @= sparse_grid_prop_sol
    sp_prop_sol[:, :, :, :] *= sp_fu1[:, None, None, :].conj()
    if is_test():
        q.json_results_append(
            f"sp_prop_sol {idx}", q.get_data_sig_arr(sp_prop_sol, rs_sig, 3), 1e-7
        )
    assert isinstance(sp_prop_sol, q.PselProp), f"sparse_solve: sp_prop_sol type={type(sp_prop_sol)}"
    return sp_prop_sol

# --------------------------------------------

@q.timer(is_timer_fork=True)
def gen_test_data():
    """
    Build a synthetic command-line argument list for testing.\n
    Uses the first entry from ``job_tag_list`` (``"test-4nt8-checker"``), generates
    a sample gauge field via the pipeline in ``qlat_scripts``, and returns an ``argv``
    list with ``--gf``, ``--out``, and all solver parameters suitable for a quick test.\n
    Returns
    -------
    argv : list of str
        Simulated command-line arguments for ``run()``.
    """
    job_tag_list = [
        "test-4nt8-checker",
        # "test-8nt16-checker",
    ]
    job_tag_traj_list = []
    for job_tag in job_tag_list:
        run_params(job_tag)
        traj_list = get_param(job_tag, "traj_list")
        for traj in traj_list:
            job_tag_traj_list.append(
                (
                    job_tag,
                    traj,
                )
            )
    fn_gf_list = []
    fn_out_list = []
    for job_tag, traj in job_tag_traj_list:
        traj_gf = traj
        run_gf(job_tag, traj_gf)
        fn_gf = get_load_path(
            f"{job_tag}/configs/ckpoint_lat.{traj_gf}",
            f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",
        )
        assert fn_gf is not None
        fn_out = get_save_path(f"{job_tag}/topo-dwf-info/traj-{traj_gf}")
        fn_gf_list.append(fn_gf)
        fn_out_list.append(fn_out)
    argv = []
    for fn in fn_gf_list:
        argv += [
            "--gf",
            fn,
        ]
    for fn in fn_out_list:
        argv += [
            "--out",
            fn,
        ]
    argv += [
        "--sparse_ratio",
        "2",
        "--num_of_rand_vol_u1",
        "2",
        "--ls",
        "16",
        "--ls_sloppy",
        "12",
        "--b_plus_c",
        "3.0",
        "--b_plus_c_sloppy",
        "3.0",
        "--maxiter_sloppy",
        "50",
        "--maxiter_exact",
        "100",
    ]
    assert isinstance(argv, list), f"gen_test_data: argv type={type(argv)}"
    for a in argv:
        assert isinstance(a, str), f"gen_test_data: argv element type={type(a)}"
    return argv

@q.timer(is_timer_fork=True)
def run_topo_measure(fn_gf, fn_out, *, params=None):
    """
    Run the topological charge measurement for a single gauge field file.\n
    Skips if the output checkpoint already exists or if the gauge field file is missing.
    Acquires a file lock to prevent concurrent execution on the same output directory.\n
    Parameters
    ----------
    fn_gf : str
        Path to the input gauge field file.
    fn_out : str
        Output directory for measurement results.
    params : dict or None
        Parameter dictionary forwarded to ``measure_topo_dwf``.  If ``"seed"`` is not
        present, it is set to ``fn_gf``.
    """
    fname = q.get_fname()
    q.check_time_limit()
    assert isinstance(fn_gf, str), f"run_topo_measure: fn_gf type={type(fn_gf)}"
    assert isinstance(fn_out, str), f"run_topo_measure: fn_out type={type(fn_out)}"
    assert params is None or isinstance(params, dict), f"run_topo_measure: params type={type(params)}"
    if q.does_file_exist_qar_sync_node(fn_out + "/checkpoint.txt"):
        q.displayln_info(
            -1, f"{fname}: WARNING: '{fn_out}' for '{fn_gf}' already exist. Skip."
        )
        return
    if not q.does_file_exist_qar_sync_node(fn_gf):
        q.displayln_info(
            -1, f"{fname}: WARNING: '{fn_gf}' does not exist. Skip this file."
        )
        return
    if not q.obtain_lock(f"{fn_out}/locks/measure-topo-dwf"):
        return
    q.json_results_append(
        f"{fname}: Start compute topo_measure out='{fn_out}' for '{fn_gf}'"
    )
    gf = q.GaugeField()
    gf.load(fn_gf)
    if "seed" not in params:
        params = dict(
            params,
            seed=fn_gf,
        )
    measure_topo_dwf(
        gf,
        info_path=fn_out,
        params=params,
    )
    q.release_lock()
    q.qremove_all_info(f"{fn_out}/locks")
    q.json_results_append(
        f"{fname}: End compute topo_measure out='{fn_out}' for '{fn_gf}'"
    )

def show_usage():
    """Print the usage / help message to stdout and exit."""
    q.displayln_info(f"Usage:{usage}")

@q.timer(is_timer_fork=True)
def run():
    """
    Main entry point for the topological charge measurement workflow.\n
    In test mode (``--test``), generates test data via ``gen_test_data()`` and uses
    the resulting synthetic arguments.  Otherwise parses ``sys.argv``.\n
    Iterates over all ``--gf`` / ``--out`` pairs and calls ``run_topo_measure`` on each.\n
    Parameters are extracted from the command line with sensible defaults:
      ``--sparse_ratio`` (32), ``--num_of_rand_vol_u1`` (2), ``--ls`` (16),
      ``--ls_sloppy`` (12), ``--b_plus_c`` (3.0), ``--b_plus_c_sloppy`` (3.0),
      ``--maxiter_sloppy`` (50), ``--maxiter_exact`` (100), ``--ama_prob`` (0.1).
    """
    q.check_time_limit()
    if is_test():
        q.displayln_info("Will now generate test data and run topo measure.")
        argv = gen_test_data()
    else:
        argv = sys.argv
    fn_gf_list = q.get_arg_list("--gf", argv=argv)
    fn_out_list = q.get_arg_list("--out", argv=argv)
    assert len(fn_gf_list) == len(fn_out_list)
    params = dict(
        sparse_ratio=int(q.get_arg("--sparse_ratio", "32", argv=argv)),
        num_of_rand_vol_u1=int(q.get_arg("--num_of_rand_vol_u1", "2", argv=argv)),
        ls=int(q.get_arg("--ls", "16", argv=argv)),
        ls_sloppy=int(q.get_arg("--ls_sloppy", "12", argv=argv)),
        b_plus_c=float(q.get_arg("--b_plus_c", "3.0", argv=argv)),
        b_plus_c_sloppy=float(q.get_arg("--b_plus_c_sloppy", "3.0", argv=argv)),
        maxiter_sloppy=int(q.get_arg("--maxiter_sloppy", "50", argv=argv)),
        maxiter_exact=int(q.get_arg("--maxiter_exact", "100", argv=argv)),
        ama_prob=float(q.get_arg("--ama_prob", "0.1", argv=argv)),
    )
    for (
        fn_gf,
        fn_out,
    ) in zip(fn_gf_list, fn_out_list):
        q.check_time_limit()
        run_topo_measure(
            fn_gf,
            fn_out,
            params=params,
        )

# --------------------------------------------

job_tag = "test-4nt8-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        4,
        4,
        4,
        8,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(20)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(3.0)

job_tag = "test-8nt16-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        8,
        8,
        8,
        16,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(10)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-16nt32-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        16,
        16,
        16,
        32,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-48nt96-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        48,
        48,
        48,
        96,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-64nt128-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        64,
        64,
        64,
        64,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

# --------------------------------------------

if __name__ == "__main__":
    q.check_time_limit()
    is_show_usage = q.get_option("--usage")
    if is_show_usage:
        show_usage()
        exit()
    qg.begin_with_gpt()
    get_all_cexpr()
    run()
    q.timer_display()
    if is_test():
        q.check_log_json(__file__, check_eps=1e-10)
    qg.end_with_gpt()
    q.displayln_info("CHECK: finished successfully.")
