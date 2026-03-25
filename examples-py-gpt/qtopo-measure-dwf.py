#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g
import numpy as np
import sys

from auto_contractor import (
    mk_fac,
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
    return ``f_tadpole_loop_sum_arr``
    ``
    f_tadpole_loop_sum_arr[rand_vol_u1_idx, tslice_idx, 0] == topo_tslice_sum
    f_tadpole_loop_sum_arr[rand_vol_u1_idx, tslice_idx, 1] == quark_condenstate_tslice_sum
    ``
    """
    q.check_time_limit()
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
    cg_f = inv.split(inv.cg({"eps": 1e-8, "maxiter": maxiter_sloppy}), mpi_split=mpi_split)
    cg = inv.split(inv.cg({"eps": 1e-8, "maxiter": maxiter_exact}), mpi_split=mpi_split)
    pc = g.qcd.fermion.preconditioner
    #
    slv_5d = inv.defect_correcting(
        inv.mixed_precision(
            inv.preconditioned(pc.eo2_ne(), cg),
            g.single, g.double,
        ),
        eps=1e-8, maxiter=100,
    )
    slv_5d_f = inv.mixed_precision(
        inv.preconditioned(pc.eo2_ne(), cg_f),
        g.single, g.double,
    )
    #
    gf = gf.copy()
    gf.unitarize()
    #
    geo = gf.geo
    total_site = geo.total_site
    #
    q.json_results_append(f"gf.plaq()", gf.plaq(), 1e-12)
    q.json_results_append(f"gf.link_trace()", gf.link_trace(), 1e-12)
    #
    gt = q.GaugeTransform(geo)
    gt.set_rand(rs.split("gt-init"), 0.5, 50)
    gt.unitarize()
    gf = gt * gf
    q.json_results_append(f"after transform gf.plaq()", gf.plaq(), 1e-12)
    q.json_results_append(f"after transform gf.link_trace()", gf.link_trace(), 1e-12)
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
        q.json_results_append(f"{sparse_ratio=} {idx=} {len(psel)=} {q.hash_sha256(psel)=}")
    #
    q.json_results_append(f"{len(psel_list)=}")
    #
    gpt_gf = qg.gpt_from_qlat(gf)
    q.json_results_append(f"g.qcd.gauge.plaquette(gpt_gf)", g.qcd.gauge.plaquette(gpt_gf), 1e-12)
    #
    gf1 = qg.qlat_from_gpt(gpt_gf)
    q.json_results_append(f"gf1.plaq()", gf1.plaq(), 1e-12)
    #
    qm = g.qcd.fermion.mobius(gpt_gf, mobius_params)
    qm_f = g.qcd.fermion.mobius(gpt_gf, mobius_params_sloppy)
    slv_qm = qm.propagator(slv_5d).grouped(n_grouped)
    slv_qm_f = qm_f.propagator(slv_5d_f).grouped(n_grouped)
    inv_qm = qg.InverterGPT(inverter=slv_qm, qtimer=q.Timer("py:slv_qm", True))
    inv_qm_f = qg.InverterGPT(inverter=slv_qm_f, qtimer=q.Timer("py:slv_qm_f", True))
    #
    rs_rand_u1 = rs.split(f"qtopo-measure-dwf(rand_u1)")
    rs_ama = rs.split(f"qtopo-measure-dwf(ama)")
    #
    f_tadpole_loop_sum_list = []
    #
    for rand_vol_u1_idx in range(num_of_rand_vol_u1):
        q.check_time_limit()
        fu1 = q.mk_rand_vol_u1(geo, 1, rs_rand_u1.split(f"{rand_vol_u1_idx}"))
        prop_src = q.Prop(geo)
        prop_src.set_unit()
        prop_src[:, :, :, :] *= fu1[:, None, None, :]
        #
        q.json_results_append(f"fu1", q.get_data_sig_arr(fu1, rs_sig, 3), 1e-12)
        q.json_results_append(f"prop_src", q.get_data_sig_arr(prop_src, rs_sig, 3), 1e-12)
        #
        def mk_path(idx):
            return f"{info_path}/scratch/rand_vol_u1_idx-{rand_vol_u1_idx}/sparse_solve_idx-{idx}"
        #
        def check(idx):
            if not q.does_file_exist_qar_sync_node(f"{mk_path(idx)}/psel.lati"):
                return False
            if not q.does_file_exist_qar_sync_node(f"{mk_path(idx)}/sp_prop_sol.lat"):
                return False
            return True
        #
        @q.timer(is_verbose=True)
        def save(idx, sp_prop_sol):
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
            return ``sp_prop_sol_il``
            Note: ``sp_prop_sol_il`` is a list ``q.PselProp`` with ``"l"`` as the ``points_dist_type``.
            """
            id_node_list_for_shuffle = q.get_id_node_list_for_shuffle()
            root_il = [
                id_node_list_for_shuffle[
                    i % len(id_node_list_for_shuffle)
                ]
                for i, idx in enumerate(idx_list)
            ]
            geo_il = [geo for idx in idx_list]
            psel_il = [psel_list[idx] for idx in idx_list]
            ssp = q.SelectedShufflePlan("g_from_l", psel_il, root_il, geo_il)
            psel_g_il = ssp.psel_recv_list
            sp_prop_sol_g_il = [
                q.PselProp(psel_g_il[i])
                for i, idx in enumerate(idx_list)
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
                q.PselProp, sp_prop_sol_g_il, is_reverse=True,
            )
            assert len(sp_prop_sol_il) == len(idx_list)
            return sp_prop_sol_il
        #
        for idx, psel in enumerate(psel_list):
            q.check_time_limit()
            q.json_results_append(f"sparse_solve: {idx+1}/{len(psel_list)}")
            if check(idx):
                continue
            sp_prop_sol = sparse_solve(idx, psel, prop_src, fu1, inv_qm_f)
            rs_ama_idx = rs_ama.split(f"{rand_vol_u1_idx} {idx}")
            r = rs_ama_idx.u_rand_gen()
            if r <= ama_prob:
                sp_prop_sol_ama = sparse_solve(idx, psel, prop_src, fu1, inv_qm)
                sp_prop_sol_ama -= sp_prop_sol
                q.json_results_append(
                    f"sp_prop_sol_ama ratio ({idx+1}/{len(psel_list)}) ({ama_prob=:.4f})",
                    np.sqrt(q.glb_sum(q.qnorm(sp_prop_sol_ama)) /
                            q.glb_sum(q.qnorm(sp_prop_sol))).item(),
                    1e-5,
                )
                sp_prop_sol_ama *= 1 / ama_prob
                sp_prop_sol += sp_prop_sol_ama
            save(idx, sp_prop_sol)
            sp_prop_sol_il = load([idx,])
            assert np.all(
                q.get_data_sig_arr(sp_prop_sol, rs_sig, 3)
                ==
                q.get_data_sig_arr(sp_prop_sol_il[0], rs_sig, 3)
            )
        #
        sp_prop_sol_list = load([idx for idx in range(len(psel_list))])
        prop_sol = q.Prop(geo)
        prop_sol.set_zero()
        for sp_prop_sol in sp_prop_sol_list:
            prop_sol @= sp_prop_sol
        #
        q.json_results_append(f"prop_sol", q.get_data_sig_arr(prop_sol, rs_sig, 3), 1e-7)
        #
        def get_prop(flavor, p_snk, p_src):
            assert flavor == "c"
            assert isinstance(p_snk, tuple) and isinstance(p_src, tuple)
            type_snk, pos_snk = p_snk
            type_src, pos_src = p_src
            assert type_snk == "point-snk"
            assert type_src == "point-snk"
            pos_snk = q.Coordinate(pos_snk)
            pos_src = q.Coordinate(pos_src)
            assert pos_snk == pos_src
            index = geo.index_from_g_coordinate(pos_snk)
            return prop_sol.get_elem_wm(index)
        #
        f_tadpole_loop = q.FieldComplexD(geo, 2)
        f_tadpole_loop.set_zero()
        #
        cexpr = get_cexpr_tadpole_loop()
        expr_names = get_expr_names(cexpr)
        chunk_list = q.get_chunk_list(
            xg_arr_full, chunk_number=max(1, q.get_q_num_mp_processes()))
        @q.timer(is_verbose=True)
        def eval(chunk_idx):
            chunk = chunk_list[chunk_idx]
            val_arr = np.zeros((len(chunk), len(expr_names),), dtype=np.complex128)
            for idx, xg in enumerate(chunk):
                pd = {
                    "x_1": ("point-snk", tuple(xg))
                }
                val_arr[idx] = eval_cexpr(cexpr, positions_dict=pd, get_prop=get_prop)
            return val_arr
        #
        val_arr_list = q.parallel_map(eval, range(len(chunk_list)))
        #
        for idx, chunk in enumerate(chunk_list):
            f_tadpole_loop.set_elems_xg(chunk, val_arr_list[idx])
        #
        q.json_results_append(f"f_tadpole_loop sig", q.get_data_sig_arr(f_tadpole_loop, rs_sig, 3), 1e-7)
        #
        f_tadpole_loop_sum = f_tadpole_loop.glb_sum_tslice()[:]
        #
        prop_sol.save_float_from_double(
            f"{info_path}/field/rand_vol_u1_idx-{rand_vol_u1_idx}/prop_sol.field",
        )
        f_tadpole_loop.save_double(
            f"{info_path}/field/rand_vol_u1_idx-{rand_vol_u1_idx}/f_tadpole_loop.field",
        )
        q.save_pickle_obj(
            f_tadpole_loop_sum, f"{info_path}/pickle/rand_vol_u1_idx-{rand_vol_u1_idx}/f_tadpole_loop_sum.pickle")
        #
        q.save_json_obj(dict(
            total=f_tadpole_loop_sum[:, 0].real.sum().item(),
            tslice_sum=f_tadpole_loop_sum[:, 0].real,
        ), f"{info_path}/info/rand_vol_u1_idx-{rand_vol_u1_idx}/topo.json",
        )
        q.save_json_obj(dict(
            total=f_tadpole_loop_sum[:, 1].real.sum().item(),
            tslice_sum=f_tadpole_loop_sum[:, 1].real,
        ), f"{info_path}/info/rand_vol_u1_idx-{rand_vol_u1_idx}/quark_condensate.json",
        )
        #
        q.json_results_append(
            f"f_tadpole_loop_sum sig",
            q.get_data_sig_arr(f_tadpole_loop_sum, rs_sig, 3),
            1e-7,
        )
        q.json_results_append(
            f"f_tadpole_loop_sum topo sum",
            f_tadpole_loop_sum[:, 0].sum(-1),
            1e-7,
        )
        #
        f_tadpole_loop_sum_list.append(f_tadpole_loop_sum)
    #
    f_tadpole_loop_sum_arr = np.array(
        f_tadpole_loop_sum_list, dtype=np.complex128)
    #
    q.save_json_obj(dict(
        total=f_tadpole_loop_sum_arr[:, :, 0].real.sum(-1),
        tslice_sum=f_tadpole_loop_sum_arr[:, :, 0].real,
    ), f"{info_path}/info/topo.json",
    )
    q.save_json_obj(dict(
        total=f_tadpole_loop_sum_arr[:, :, 1].real.sum(-1),
        tslice_sum=f_tadpole_loop_sum_arr[:, :, 1].real,
    ), f"{info_path}/info/quark_condensate.json",
    )
    q.save_pickle_obj(
        f_tadpole_loop_sum_arr,
        f"{info_path}/pickle/f_tadpole_loop_sum_arr.pickle",
    )
    #
    q.json_results_append(
        f"f_tadpole_loop_sum_arr sig",
        q.get_data_sig_arr(f_tadpole_loop_sum_arr, rs_sig, 3),
        1e-7,
    )
    q.json_results_append(
        f"f_tadpole_loop_sum_arr topo sum",
        f_tadpole_loop_sum_arr[:, :, 0].sum(-1),
        1e-7,
    )
    #
    topo_tslice_arr = f_tadpole_loop_sum_arr[:, :, 0].copy()
    #
    q.qtouch_info(f"{info_path}/checkpoint.txt")
    q.qremove_all_info(f"{info_path}/scratch")
    #
    return f_tadpole_loop_sum_arr

# ------------------------------

@q.timer
def get_cexpr_tadpole_loop():
    fn_base = "cache/auto_contract_cexpr/get_cexpr_tadpole_loop"
    def calc_cexpr():
        diagram_type_dict = dict()
        diagram_type_dict[((('x_1', 'x_1'), 1))] = "Type1"
        exprs = [
            mk_scalar5("c", "c", "x_1") + f"qbar g5 q",
            mk_scalar("c", "c", "x_1") + f"qbar q",
        ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit=True, diagram_type_dict=diagram_type_dict)
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, fn_base, is_cython=is_cython)

@q.timer(is_verbose=True)
def get_all_cexpr():
    benchmark_eval_cexpr(get_cexpr_tadpole_loop())

@q.timer(is_verbose=True)
def mk_sparse_grid(xg_arr, sparse_ratio, idx):
    if sparse_ratio == 1:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        sel_arr &= (xg_arr[:, 0] + xg_arr[:, 1] +
                    xg_arr[:, 2] + xg_arr[:, 3]) % 1 == 0
        return sel_arr
    elif sparse_ratio == 2:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        sel_arr &= (xg_arr[:, 0] + xg_arr[:, 1] +
                    xg_arr[:, 2] + xg_arr[:, 3] + idx0) % 2 == 0
        sel2_arr = sel_arr.copy()
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
    elif sparse_ratio == 4:
        assert 0 <= idx < sparse_ratio
        sel_arr = True
        idx0 = idx % 2
        idx1 = (idx // 2) % 2
        sel_arr &= (xg_arr[:, 0] + xg_arr[:, 1] + idx0) % 2 == 0
        sel_arr &= (xg_arr[:, 2] + xg_arr[:, 3] + idx1) % 2 == 0
        sel2_arr = sel_arr.copy()
        sel2_arr &= (xg_arr[:, 0] + xg_arr[:, 1] +
                     xg_arr[:, 2] + xg_arr[:, 3] + idx0 + idx1) % 2 == 0
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        sel2_arr &= (xg_arr[:, 0] + xg_arr[:, 1] +
                     xg_arr[:, 2] + xg_arr[:, 3] + idx0 + idx1) % 2 == 0
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        sel_arr &= ((xg_arr[:, 0] + idx0) // 2 + (xg_arr[:, 1] + idx1) // 2 +
                    (xg_arr[:, 2] + idx2) // 2 + (xg_arr[:, 3] + idx3) // 2 + idx4) % 2 == 0
        sel2_arr = sel_arr.copy()
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        sel_arr &= ((xg_arr[:, 0] + idx0) // 2 +
                    (xg_arr[:, 1] + idx1) // 2 + idx4) % 2 == 0
        sel_arr &= ((xg_arr[:, 2] + idx2) // 2 +
                    (xg_arr[:, 3] + idx3) // 2 + idx5) % 2 == 0
        sel2_arr = sel_arr.copy()
        sel2_arr &= ((xg_arr[:, 0] + idx0) // 2 + (xg_arr[:, 1] + idx1) // 2 +
                     (xg_arr[:, 2] + idx2) // 2 + (xg_arr[:, 3] + idx3) // 2 + idx4 + idx5) % 2 == 0
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        sel_arr &= ((xg_arr[:, 0] + idx0) // 2 +
                    (xg_arr[:, 1] + idx1) // 2 + idx4) % 2 == 0
        sel_arr &= ((xg_arr[:, 2] + idx2) // 2 +
                    (xg_arr[:, 3] + idx3) // 2 + idx5) % 2 == 0
        sel_arr &= ((xg_arr[:, 0] + idx0) // 2 +
                    (xg_arr[:, 2] + idx2) // 2 + idx6) % 2 == 0
        sel2_arr = sel_arr.copy()
        sel2_arr &= ((xg_arr[:, 1] + idx1) // 2 -
                     (xg_arr[:, 2] + idx2) // 2 + idx4 - idx6) % 2 == 0
        sel2_arr &= ((xg_arr[:, 0] + idx0) // 2 -
                     (xg_arr[:, 3] + idx3) // 2 + idx6 - idx5) % 2 == 0
        sel2_arr &= ((xg_arr[:, 0] + idx0) // 2 + (xg_arr[:, 1] + idx1) // 2 +
                     (xg_arr[:, 2] + idx2) // 2 + (xg_arr[:, 3] + idx3) // 2 + idx4 + idx5) % 2 == 0
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        sel_arr &= ((xg_arr[:, 0] + idx0) // 3 + (xg_arr[:, 1] + idx1) // 3 +
                    (xg_arr[:, 2] + idx2) // 3 + (xg_arr[:, 3] + idx3) // 3 + idx4) % 2 == 0
        sel2_arr = sel_arr.copy()
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
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
        sel_arr &= ((xg_arr[:, 0] + idx0) // 4 + (xg_arr[:, 1] + idx1) // 4 +
                    (xg_arr[:, 2] + idx2) // 4 + (xg_arr[:, 3] + idx3) // 4 + idx4) % 2 == 0
        sel2_arr = sel_arr.copy()
        assert np.all(sel2_arr == sel_arr)
        return sel_arr
    else:
        raise Exception(f"{sparse_ratio=}")

@q.timer(is_verbose=True)
def sparse_solve(idx, psel, prop_src, fu1, inverter):
    """
    return ``sp_prop_sol``
    Complex conjugate of the initial random phase is already multiplied.
    """
    geo = prop_src.geo
    sp_prop_src = q.PselProp(psel)
    sp_prop_src.set_zero()
    sp_prop_src @= prop_src
    if is_test():
        q.json_results_append(f"sp_prop_src {idx}", q.get_data_sig_arr(sp_prop_src, rs_sig, 3), 1e-12)
    sparse_grid_prop_src = q.Prop(geo)
    sparse_grid_prop_src.set_zero()
    sparse_grid_prop_src @= sp_prop_src
    sparse_grid_prop_sol = inverter * sparse_grid_prop_src
    sp_prop_sol = q.PselProp(psel)
    sp_prop_sol @= sparse_grid_prop_sol
    sp_fu1 = q.SelectedPointsComplexD(psel, 1)
    sp_fu1 @= fu1
    sp_prop_sol[:, :, :, :] *= sp_fu1[:, None, None, :].conj()
    if is_test():
        q.json_results_append(f"sp_prop_sol {idx}", q.get_data_sig_arr(sp_prop_sol, rs_sig, 3), 1e-7)
    return sp_prop_sol

# --------------------------------------------

@q.timer(is_timer_fork=True)
def gen_test_data():
    job_tag_list = [
            "test-4nt8-checker",
            # "test-8nt16-checker",
            ]
    job_tag_traj_list = []
    for job_tag in job_tag_list:
        run_params(job_tag)
        traj_list = get_param(job_tag, "traj_list")
        for traj in traj_list:
            job_tag_traj_list.append((job_tag, traj,))
    fn_gf_list = []
    fn_out_list = []
    for job_tag, traj in job_tag_traj_list:
        traj_gf = traj
        get_gf = run_gf(job_tag, traj_gf)
        fn_gf = get_load_path(f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}")
        assert fn_gf is not None
        fn_out = get_save_path(f"{job_tag}/topo-dwf-info/traj-{traj_gf}")
        fn_gf_list.append(fn_gf)
        fn_out_list.append(fn_out)
    argv = []
    for fn in fn_gf_list:
        argv += [ "--gf", fn, ]
    for fn in fn_out_list:
        argv += [ "--out", fn, ]
    argv += [
        "--sparse_ratio", "2",
        "--num_of_rand_vol_u1", "2",
        "--ls", "16",
        "--ls_sloppy", "12",
        "--b_plus_c", "3.0",
        "--b_plus_c_sloppy", "3.0",
    ]
    return argv

@q.timer(is_timer_fork=True)
def run_topo_measure(fn_gf, fn_out, *, params=None):
    fname = q.get_fname()
    q.check_time_limit()
    if q.does_file_exist_qar_sync_node(fn_out + "/checkpoint.txt"):
        q.displayln_info(-1, f"{fname}: WARNING: '{fn_out}' for '{fn_gf}' already exist. Skip.")
        return
    if not q.does_file_exist_qar_sync_node(fn_gf):
        q.displayln_info(-1, f"{fname}: WARNING: '{fn_gf}' does not exist. Skip this file.")
        return
    if not q.obtain_lock(f"{fn_out}/locks/measure-topo-dwf"):
        return
    q.json_results_append(f"{fname}: Start compute topo_measure out='{fn_out}' for '{fn_gf}'")
    gf = q.GaugeField()
    gf.load(fn_gf)
    if "seed" not in params:
        params = dict(
            params,
            seed=fn_gf,
        )
    f_tadpole_loop_sum_arr = measure_topo_dwf(
        gf,
        info_path=fn_out,
        params=params,
    )
    q.release_lock()
    q.qremove_all_info(f"{fn_out}/locks")
    q.json_results_append(f"{fname}: End compute topo_measure out='{fn_out}' for '{fn_gf}'")

def show_usage():
    q.displayln_info(f"Usage:{usage}")

@q.timer(is_timer_fork=True)
def run():
    q.check_time_limit()
    if is_test():
        q.displayln_info(f"Will now generate test data and run topo measure.")
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
    for fn_gf, fn_out, in zip(fn_gf_list, fn_out_list):
        q.check_time_limit()
        run_topo_measure(
            fn_gf,
            fn_out,
            params=params,
        )

# --------------------------------------------

job_tag = "test-4nt8-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(20)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(3.0)

job_tag = "test-8nt16-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 8, 8, 8, 16, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(10)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-16nt32-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 16, 16, 16, 32, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-48nt96-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 48, 48, 48, 96, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-64nt128-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 64, 64, 64, 64, ])
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
    q.displayln_info(f"CHECK: finished successfully.")