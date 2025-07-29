#!/usr/bin/env python3

import qlat as q
import numpy as np
import sys

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin_with_mpi(size_node_list)

total_site_list = [
        q.Coordinate([ 4, 4, 4, 4, ]),
        q.Coordinate([ 6, 6, 6, 6, ]),
        q.Coordinate([ 8, 8, 8, 8, ]),
        ]

multiplicity_list = [
        1, 2, 3,
        ]

def get_f_list_sig(f_list, rs, n):
    sig = np.zeros(n, dtype=np.complex128)
    for idx, sp in enumerate(f_list):
        sig += q.get_data_sig_arr(sp[:], rs.split(f"sig {idx}"), n)
    sig = q.glb_sum(sig)
    return sig

@q.timer
def selected_shuffle_r_from_l(total_site, multiplicity, seed):
    fname = q.get_fname()
    q.json_results_append(f"{fname}: {total_site} {multiplicity} {seed}")
    rs = q.RngState(f"seed {fname} {seed}")
    #
    n_points = total_site.volume() // 16
    q.json_results_append(f"n_points={n_points}")
    #
    geo = q.Geometry(total_site)
    psel = q.PointsSelection()
    psel.set_rand(total_site, n_points, rs.split("psel"))
    q.json_results_append(f"hash(psel)={q.hash_sha256(psel)}")
    fsel = q.FieldSelection(psel)
    psel_l = q.PointsSelection(fsel)
    #
    ssp = q.SelectedShufflePlan("r_from_l", psel_l, geo, rs.split("ssp"))
    assert psel_l is ssp.psel_src_list[0]
    psel_s = ssp.psel_dst_list[0]
    #
    psel_s1 = q.PointsSelection(psel_l, ssp)
    assert psel_s1 == psel_s
    #
    q.json_results_append(f"hash(psel_s)={q.hash_sha256(psel_s)}")
    q.displayln_info(f"len(psel)={len(psel)} ; psel={psel}")
    psel_str = f"len(psel_s)={len(psel_s)} ; psel_s={psel_s}"
    psel_str_list = q.get_comm().allgather(psel_str)
    for id_node, psel_str in enumerate(psel_str_list):
        q.displayln_info(f"id_node={id_node} ; {psel_str}")
    #
    sp_l = q.SelectedPointsComplexD(psel_l, multiplicity)
    assert sp_l.psel is ssp.psel_src_list[0]
    sp_l.set_rand(rs.split("sp_l"))
    sig_l = q.get_data_sig_arr(sp_l, rs.split("sig"), 2)
    q.json_results_append(f"get_data_sig_arr(sp_l,rs,2)", sig_l, 1e-12)
    sp_s = q.SelectedPointsComplexD(sp_l, ssp)
    assert sp_s.psel is ssp.psel_dst_list[0]
    sig_s = q.get_data_sig_arr(sp_s, rs.split("sig"), 2)
    q.json_results_append(f"get_data_sig_arr(sp_s,rs,2)", sig_s, 1e-12)
    psel_l1 = q.PointsSelection(psel_s, ssp, True)
    assert psel_l1 == psel_l
    sp_l1 = q.SelectedPointsComplexD(sp_s, ssp, True)
    assert sp_l1.psel is ssp.psel_src_list[0]
    sig_l1 = q.get_data_sig_arr(sp_l1, rs.split("sig"), 2)
    q.json_results_append(f"get_data_sig_arr(sp_l1,rs,2)", sig_l1, 1e-12)
    assert np.all(sig_l1 == sig_l)
    #
    sp_lc = q.SelectedPointsChar()
    sp_l1.swap_cast(sp_lc)
    sp_sc = q.SelectedPointsChar(sp_lc, ssp)
    sp_l1.swap_cast(sp_sc)
    sp_s1 = sp_l1
    sig_s1 = q.get_data_sig_arr(sp_s1, rs.split("sig"), 2)
    q.json_results_append(f"get_data_sig_arr(sp_s1,rs,2)", sig_s1, 1e-12)
    assert np.all(sig_s1 == sig_s)
    #
    psel_l_list = [ psel_l.copy() for i in range(3) ]
    geo_list = [ geo.copy() for i in range(3) ]
    psel_l_list_hash = q.hash_sha256(psel_l_list)
    q.json_results_append(f"hash(psel_l_list)={psel_l_list_hash}")
    #
    ssp = q.SelectedShufflePlan("r_from_l", psel_l_list, geo_list, rs.split("ssp"))
    assert psel_l_list is ssp.psel_src_list
    assert psel_l_list_hash == q.hash_sha256(ssp.psel_src_list)
    q.json_results_append(f"hash(ssp.psel_dst_list)={q.hash_sha256(ssp.psel_dst_list)}")
    #
    sp_l_list = []
    for idx, psel in enumerate(psel_l_list):
        sp = q.SelectedPointsComplexD(psel, multiplicity)
        sp.set_rand(rs.split(f"sp_l {idx}"))
        sp_l_list.append(sp)
    q.json_results_append(f"sig sp_l_list", get_f_list_sig(sp_l_list, rs, 3), 1e-12)
    #
    sp_s_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_l_list)
    q.json_results_append(f"sig sp_s_list", get_f_list_sig(sp_s_list, rs, 3), 1e-12)
    sp_ss_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_s_list, is_reverse=True)
    assert np.all(get_f_list_sig(sp_ss_list, rs, 3) == get_f_list_sig(sp_l_list, rs, 3))

@q.timer
def selected_shuffle_t_slice_from_l(total_site, multiplicity, seed):
    fname = q.get_fname()
    q.json_results_append(f"{fname}: {total_site} {multiplicity} {seed}")
    rs = q.RngState(f"seed {fname} {seed}")
    geo = q.Geometry(total_site)
    #
    n_points = total_site.volume() // 16
    q.json_results_append(f"n_points={n_points}")
    #
    geo = q.Geometry(total_site)
    psel = q.PointsSelection()
    psel.set_rand(total_site, n_points, rs.split("psel"))
    q.json_results_append(f"hash(psel)={q.hash_sha256(psel)}")
    fsel = q.FieldSelection(psel)
    psel_l = q.PointsSelection(fsel)
    #
    psel_l_list = [ psel_l.copy() for i in range(3) ]
    geo_list = [ geo.copy() for i in range(3) ]
    psel_l_list_hash = q.hash_sha256(psel_l_list)
    #
    q.json_results_append(f"hash(psel_l_list)={psel_l_list_hash}")
    ssp = q.SelectedShufflePlan("t_slice_from_l", psel_l_list, geo_list)
    assert psel_l_list is ssp.psel_src_list
    assert psel_l_list_hash == q.hash_sha256(ssp.psel_src_list)
    q.json_results_append(f"hash(ssp.psel_dst_list)={q.hash_sha256(ssp.psel_dst_list)}")
    #
    sp_l_list = []
    for idx, psel in enumerate(psel_l_list):
        sp = q.SelectedPointsComplexD(psel, multiplicity)
        sp.set_rand(rs.split(f"sp_l {idx}"))
        sp_l_list.append(sp)
    q.json_results_append(f"sig sp_l_list", get_f_list_sig(sp_l_list, rs, 3), 1e-12)
    #
    sp_s_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_l_list)
    q.json_results_append(f"sig sp_s_list", get_f_list_sig(sp_s_list, rs, 3), 1e-12)
    sp_ss_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_s_list, is_reverse=True)
    assert np.all(get_f_list_sig(sp_ss_list, rs, 3) == get_f_list_sig(sp_l_list, rs, 3))
    #
    ssp = q.SelectedShufflePlan("dist_t_slice_from_l", psel_l, geo, len(psel_l_list))
    assert psel_l is ssp.psel_src_list[0]
    q.json_results_append(f"hash(ssp.psel_dst_list)={q.hash_sha256(ssp.psel_dst_list)}")
    #
    sp_l = q.SelectedPointsComplexD(psel_l, multiplicity)
    assert sp_l.psel is ssp.psel_src_list[0]
    sp_l.set_rand(rs.split("sp_l"))
    sig_l = q.get_data_sig_arr(sp_l, rs.split("sig"), 2)
    q.json_results_append(f"get_data_sig_arr(sp_l,rs,2)", sig_l, 1e-12)
    #
    sp_s_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, [ sp_l, ])
    q.json_results_append(f"sig sp_s_list", get_f_list_sig(sp_s_list, rs, 3), 1e-12)
    sp_ss_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_s_list, is_reverse=True)
    assert np.all(get_f_list_sig(sp_ss_list, rs, 3) == get_f_list_sig([ sp_l, ], rs, 3))

@q.timer
def selected_shuffle_t_slice_from_f(total_site, multiplicity, seed):
    fname = q.get_fname()
    q.json_results_append(f"{fname}: {total_site} {multiplicity} {seed}")
    rs = q.RngState(f"seed {fname} {seed}")
    geo = q.Geometry(total_site)
    psel = q.PointsSelection(geo)
    q.json_results_append(f"len(psel)={len(psel)}")
    q.json_results_append(f"hash_sha256(psel)={q.hash_sha256(psel)}")
    #
    t_slice_node_site = total_site.copy()
    t_slice_node_site[3] = 1
    t_slice_size_node = q.Coordinate([ 1, 1, 1, total_site[3], ])
    #
    num_field = 3
    #
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split(f"gf"), 0.5, 2)
    gf_sig = q.get_data_sig_arr(gf, rs, 3)
    q.json_results_append(f"gf sig", gf_sig, 1e-10)
    #
    psel_list = []
    geo_list = []
    vec_list = []
    for i in range(num_field):
        vec = q.FermionField4d(geo)
        vec.set_rand(rs.split(f"vec {i}"))
        psel_list.append(psel.copy())
        geo_list.append(geo.copy())
        vec_list.append(vec)
    q.json_results_append(f"hash_sha256(psel_list)={q.hash_sha256(psel_list)}")
    vec_sig = get_f_list_sig(vec_list, rs, 3)
    q.json_results_append(f"sig vec_list", vec_sig, 1e-12)
    #
    spc = q.SelectedPointsChar(psel)
    geo = q.Geometry()
    gf.swap_sp_cast(spc, geo)
    #
    gf = q.GaugeField()
    gf.swap_sp_cast(spc, geo)
    assert np.all(gf_sig == q.get_data_sig_arr(gf, rs, 3))
    #
    ssp1 = q.SelectedShufflePlan("dist_t_slice_from_l", psel, geo, num_field)
    assert ssp1.psel_src_list[0] is psel
    q.json_results_append(f"hash_sha256(ssp1.psel_dst_list)={q.hash_sha256(ssp1.psel_dst_list)}")
    t_list = [ psel[0][3].item() for psel in ssp1.psel_dst_list ]
    s_geo_list = [ q.Geometry(q.Coordinate([ 0, 0, 0, t, ]), t_slice_size_node, t_slice_node_site) for t in t_list ]
    all_t_list_list = q.get_comm().allgather(t_list)
    q.json_results_append(f"all_t_list_list={all_t_list_list}")
    #
    gf.swap_sp_cast(spc, geo)
    s_spc_list = ssp1.shuffle_list([ spc, ])
    s_gf_list = []
    for s_spc, s_geo in zip(s_spc_list, s_geo_list):
        s_gf = q.GaugeField()
        s_gf.swap_sp_cast(s_spc, s_geo)
        s_gf_list.append(s_gf)
    q.json_results_append(f"sig s_gf_list", get_f_list_sig(s_gf_list, rs, 3), 1e-10)
    for s_gf, s_spc, s_geo in zip(s_gf_list, s_spc_list, s_geo_list):
        s_gf.swap_sp_cast(s_spc, s_geo)
    [ spc, ] = ssp1.shuffle_list(s_spc_list, is_reverse=True)
    gf.swap_sp_cast(spc, geo)
    assert np.all(gf_sig == q.get_data_sig_arr(gf, rs, 3))
    #
    ssp2 = q.SelectedShufflePlan("t_slice_from_l", psel_list, geo_list)
    assert ssp2.psel_src_list is psel_list
    q.json_results_append(f"hash_sha256(ssp2.psel_dst_list)={q.hash_sha256(ssp2.psel_dst_list)}")
    vt_list = [ psel[0][3].item() for psel in ssp2.psel_dst_list ]
    vs_geo_list = [ q.Geometry(q.Coordinate([ 0, 0, 0, t, ]), t_slice_size_node, t_slice_node_site) for t in vt_list ]
    all_vt_list_list = q.get_comm().allgather(vt_list)
    q.json_results_append(f"all_vt_list_list={all_vt_list_list}")
    #
    vspc_list = []
    vgeo_list = []
    for i in range(num_field):
        vspc = q.SelectedPointsChar(psel)
        vgeo = q.Geometry()
        vec_list[i].swap_sp_cast(vspc, vgeo)
        vspc_list.append(vspc)
        vgeo_list.append(vgeo)
    s_vspc_list = ssp2.shuffle_list(vspc_list)
    s_vec_list = []
    for s_vspc, vs_geo in zip(s_vspc_list, vs_geo_list):
        s_vec = q.FermionField4d()
        s_vec.swap_sp_cast(s_vspc, vs_geo)
        s_vec_list.append(s_vec)
    q.json_results_append(f"sig s_vec_list", get_f_list_sig(s_vec_list, rs, 3), 1e-10)
    for s_vec, s_vspc, vs_geo in zip(s_vec_list, s_vspc_list, vs_geo_list):
        s_vec.swap_sp_cast(s_vspc, vs_geo)
    vspc_list = ssp2.shuffle_list(s_vspc_list, is_reverse=True)
    assert len(vspc_list) == num_field
    assert len(vgeo_list) == num_field
    for i in range(num_field):
        vec_list[i].swap_sp_cast(vspc_list[i], vgeo_list[i])
    assert np.all(vec_sig == get_f_list_sig(vec_list, rs, 3))

@q.timer
def test_prop_spatial_smear(total_site, multiplicity, seed):
    fname = q.get_fname()
    q.json_results_append(f"{fname}: {total_site} {multiplicity} {seed}")
    rs = q.RngState(f"seed {fname} {seed}")
    geo = q.Geometry(total_site)
    #
    coef = 0.9375
    step = 20
    #
    num_field = multiplicity
    #
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split(f"gf"), 0.5, 2)
    gf_sig = q.get_data_sig_arr(gf, rs, 3)
    q.json_results_append(f"gf sig", gf_sig, 1e-10)
    #
    ff_list = []
    for id_field in range(num_field):
        ff = q.FermionField4d(geo)
        ff.set_rand(rs.split(f"ff_list[{id_field}]"))
        ff_list.append(ff)
    ff_list_sig = get_f_list_sig(ff_list, rs, 3)
    q.json_results_append(f"sig ff_list", ff_list_sig, 1e-12)
    #
    s_ff_list = prop_spatial_smear(ff_list, gf, coef, step)
    s_ff_list_sig = get_f_list_sig(s_ff_list, rs, 3)
    q.json_results_append(f"sig s_ff_list", s_ff_list_sig, 1e-12)
    #
    mom = q.CoordinateD([ 0.1, -0.2, 0.3, 0.5, ])
    #
    ms_ff_list = prop_spatial_smear(ff_list, gf, coef, step, mom)
    ms_ff_list_sig = get_f_list_sig(ms_ff_list, rs, 3)
    q.json_results_append(f"sig ms_ff_list", ms_ff_list_sig, 1e-12)

@q.timer
def prop_spatial_smear(ff_list, gf, coef, step, mom=None):
    """
    Perform spatial smear for `ff_list`, a list of `FermionField4d`.
    Return new `ff_list` after smearing.
    Original `ff_list` should not be modified.
    #
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    gf = get_gf_ape()
    #
    coef = get_param(job_tag, "prop_smear_coef")
    step = get_param(job_tag, "prop_smear_step")
    #
    job_tag = "64I"
    set_param(job_tag, "gf_ape_smear_coef")(0.5)
    set_param(job_tag, "gf_ape_smear_step")(30)
    set_param(job_tag, "prop_smear_coef")(0.9375)
    set_param(job_tag, "prop_smear_step")(54)
    """
    fname = q.get_fname()
    #
    assert isinstance(ff_list, list)
    num_field = len(ff_list)
    for i in range(num_field):
        assert isinstance(ff_list[i], q.FermionField4d)
    assert isinstance(gf, q.GaugeField)
    assert isinstance(coef, float)
    assert isinstance(step, int)
    assert step >= 0
    if mom is None:
        mom = q.CoordinateD()
    assert isinstance(mom, q.CoordinateD)
    if step == 0:
        return
    #
    gf = gf.copy()
    geo = gf.geo
    total_site = geo.total_site
    psel = q.PointsSelection(geo)
    psel_list = [ psel.copy() for i in range(num_field) ]
    geo_list = [ geo.copy() for i in range(num_field) ]
    #
    ssp1 = q.SelectedShufflePlan("dist_t_slice_from_l", psel, geo, num_field)
    ssp2 = q.SelectedShufflePlan("t_slice_from_l", psel_list, geo_list)
    #
    s_gf_list = ssp1.shuffle_sp_list(q.GaugeField, [ gf, ])
    s_ff_list = ssp2.shuffle_sp_list(q.FermionField4d, ff_list)
    #
    id_node_set = set()
    s_ff_mask_arr = np.zeros(len(s_ff_list), dtype=np.int8)
    for s_gf in s_gf_list:
        id_node = s_gf.geo.id_node
        assert id_node not in id_node_set
        id_node_set.add(id_node)
        sub_s_ff_list = []
        for idx, s_ff in enumerate(s_ff_list):
            if id_node == s_ff.geo.id_node:
                assert s_ff_mask_arr[idx] == 0
                s_ff_mask_arr[idx] = 1
                sub_s_ff_list.append(s_ff)
        q.prop_spatial_smear_no_comm(sub_s_ff_list, s_gf, coef, step, mom)
    assert np.all(s_ff_mask_arr == 1)
    #
    ss_ff_list = ssp2.shuffle_sp_list(q.FermionField4d, s_ff_list, is_reverse=True)
    return ss_ff_list

for total_site in total_site_list:
    for multiplicity in multiplicity_list:
        for seed in range(1):
            selected_shuffle_r_from_l(total_site, multiplicity, seed)
            selected_shuffle_t_slice_from_l(total_site, multiplicity, seed)
            selected_shuffle_t_slice_from_f(total_site, multiplicity, seed)
            test_prop_spatial_smear(total_site, multiplicity, seed)

q.timer_display()
q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
