#!/usr/bin/env python3

import qlat as q
import numpy as np

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
    q.Coordinate(
        [
            4,
            4,
            4,
            4,
        ]
    ),
    q.Coordinate(
        [
            6,
            6,
            6,
            6,
        ]
    ),
    q.Coordinate(
        [
            8,
            8,
            8,
            8,
        ]
    ),
]

multiplicity_list = [
    1,
    2,
    3,
]

def get_f_list_sig(f_list, rs, n):
    sig = np.zeros(n, dtype=np.complex128)
    for idx, sp in enumerate(f_list):
        sig += q.get_data_sig_arr(sp[:], rs.split(f"sig {idx}"), n)
    sig = q.glb_sum(sig)
    return sig

@q.timer
def selected_shuffle_l_from_g(total_site, multiplicity, seed):
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
    #
    root = 0
    ssp = q.SelectedShufflePlan("l_from_g", psel, root)
    if q.get_id_node() == root:
        assert psel is ssp.psel_send_list[0]
    else:
        assert ssp.psel_send_list[0] is None
    #
    psel_l = ssp.psel_recv_list[0]
    q.json_results_append(f"hash(psel_l)={q.hash_sha256(psel_l)}")
    #
    psel_l1 = q.PointsSelection(psel, ssp)
    assert psel_l == psel_l1
    #
    psel1 = q.PointsSelection(psel_l, ssp, True)
    if q.get_id_node() == root:
        assert psel == psel1
    else:
        assert psel1.n_points == 0
        assert psel1.points_dist_type == "g"
    q.json_results_append(f"hash(psel1)={q.hash_sha256(psel1)}")
    #
    ssp_r = q.SelectedShufflePlan("g_from_l", psel_l, root, geo)
    assert psel_l is ssp_r.psel_send_list[0]
    psel_s = ssp_r.psel_recv_list[0]
    q.json_results_append(f"hash(psel_s)={q.hash_sha256(psel_s)}")
    #
    if q.get_id_node() == root:
        assert sorted(psel_s.xg_arr.tolist()) == sorted(psel1.xg_arr.tolist())
        assert psel_s.xg_arr[:, ::-1].tolist() == sorted(psel1.xg_arr[:, ::-1].tolist())
    else:
        assert psel_s.n_points == 0
        assert psel_s.points_dist_type == "g"
    #
    psel_s1 = q.PointsSelection(psel_l, ssp_r)
    assert psel_s1 == psel_s
    q.json_results_append(f"hash(psel_s1)={q.hash_sha256(psel_s1)}")
    assert q.hash_sha256(psel_s1) == q.hash_sha256(psel_s)
    #
    q.displayln_info(f"len(psel)={len(psel)} ; psel={psel}")
    q.displayln_info(f"len(psel)={len(psel_s)} ; psel={psel_s}")
    psel_str = f"len(psel_l)={len(psel_l)} ; psel_l={psel_l}"
    psel_str_list = q.get_comm().allgather(psel_str)
    for id_node, psel_str in enumerate(psel_str_list):
        q.displayln_info(f"id_node={id_node} ; {psel_str}")
    #
    sp = q.SelectedPointsComplexD(psel, multiplicity)
    if q.get_id_node() == root:
        assert sp.psel is ssp.psel_send_list[0]
    else:
        assert ssp.psel_send_list[0] is None
    sp.set_rand(rs.split("sp"))
    sig = q.get_data_sig_arr(sp, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp,rs,2)", sig, 1e-12)
    sp_l = q.SelectedPointsComplexD(sp, ssp)
    assert sp_l.psel is ssp.psel_recv_list[0]
    sig_l = q.get_data_sig_arr(sp_l, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp_l,rs,2)", sig_l, 1e-12)
    sp1 = q.SelectedPointsComplexD(sp_l, ssp, True)
    assert sp1.psel is ssp.psel_send_list[0]
    # sig1 = q.get_data_sig_arr(sp1, rs.split("sig"), 2)
    if q.get_id_node() == root:
        assert np.all(sp1[:] == sp[:])
    else:
        assert sp1.n_points == 0
        assert sp1.points_dist_type == "g"
        print("sp1[:]", sp1[:])
    q.sync_node("sp1 sig1")
    sp1.bcast(root)
    q.sync_node("sp1 bcast")
    sig1 = q.get_data_sig_arr(sp1, rs.split("sig"), 2)
    q.sync_node("sp1 sig1")
    q.json_results_append("get_data_sig_arr(sp1,rs,2)", sig1, 1e-12)
    assert np.all(sig1 == sig)
    #
    sp_l = q.SelectedPointsComplexD(psel_l, multiplicity)
    assert sp_l.psel is ssp_r.psel_send_list[0]
    sp_l.set_rand(rs.split("sp"))
    sig_l = q.get_data_sig_arr(sp_l, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp_l,rs,2)", sig_l, 1e-12)
    sp2 = q.SelectedPointsComplexD(sp_l, ssp_r)
    assert sp2.psel is ssp_r.psel_recv_list[0]
    sig2 = q.get_data_sig_arr(sp2, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp2,rs,2)", sig2, 1e-12)
    if q.get_id_node() == root:
        sp3 = q.SelectedPointsComplexD(sp2.psel, sp2.multiplicity)
        sp3 @= sp
        assert np.all(sp2[:] == sp3[:])
    else:
        assert sp2.n_points == 0
        assert sp2.points_dist_type == "g"
    sp_l1 = q.SelectedPointsComplexD(sp2, ssp_r, True)
    sig_l1 = q.get_data_sig_arr(sp_l1, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp_l1,rs,2)", sig_l1, 1e-12)
    assert np.all(sig_l1 == sig_l)
    #
    sp_l1 = sp_l.copy()
    #
    sp_lc = q.SelectedPointsChar()
    sp_l1.swap_cast(sp_lc)
    sp_sc = q.SelectedPointsChar(sp_lc, ssp_r)
    sp_l1.swap_cast(sp_sc)
    sp_s1 = sp_l1
    sig_s1 = q.get_data_sig_arr(sp_s1, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp_s1,rs,2)", sig_s1, 1e-12)
    assert np.all(sig_s1 == sig2)
    #
    psel_list = [psel.copy() for i in range(3)]
    geo_list = [geo.copy() for i in range(3)]
    root_list = [i % q.get_num_node() for i in range(3)]
    psel_list_hash = q.hash_sha256(psel_list)
    q.json_results_append(f"hash(psel_list)={psel_list_hash}")
    #
    ssp = q.SelectedShufflePlan("l_from_g", psel_list, root_list)
    ssp_r = q.SelectedShufflePlan("g_from_l", ssp.psel_recv_list, root_list, geo_list)
    for p1, p2 in zip(psel_list, ssp.psel_send_list):
        assert p1 is p2 or p2 is None
    q.json_results_append(
        f"hash(ssp.psel_recv_list)={q.hash_sha256(ssp.psel_recv_list)}"
    )
    #
    sp_list = []
    for idx, psel_idx in enumerate(psel_list):
        sp = q.SelectedPointsComplexD(psel_idx, multiplicity)
        sp.set_rand(rs.split(f"sp {idx}"))
        sp_list.append(sp)
    sig_list = get_f_list_sig(sp_list, rs, 3)
    q.json_results_append("sig sp_list", sig_list, 1e-12)
    #
    sp_l_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_list)
    sig_l_list = get_f_list_sig(sp_l_list, rs, 3)
    q.json_results_append("sig sp_l_list", sig_l_list, 1e-12)
    sp1_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_l_list, is_reverse=True)
    sig1_list = get_f_list_sig(sp1_list, rs, 3)
    q.json_results_append("sig sp1_list", sig1_list, 1e-12)
    for sp1, sp, root in zip(sp1_list, sp_list, root_list):
        if q.get_id_node() == root:
            assert np.all(sp1[:] == sp[:])
        else:
            assert sp1.n_points == 0
            assert sp1.points_dist_type == "g"
        sp1.bcast(root)
    sig1p_list = get_f_list_sig(sp1_list, rs, 3)
    q.json_results_append("sig sp1_list prime", sig1p_list, 1e-12)
    assert np.all(sig1p_list == sig_list)
    #
    sp2_list = ssp_r.shuffle_sp_list(q.SelectedPointsComplexD, sp_l_list)
    sig2_list = get_f_list_sig(sp2_list, rs, 3)
    q.json_results_append("sig sp2_list", sig2_list, 1e-12)
    #
    for sp1, sp, root in zip(sp2_list, sp_list, root_list):
        if q.get_id_node() == root:
            assert sp1.n_points == sp.n_points
            assert sp1.points_dist_type == sp.points_dist_type
        else:
            assert sp1.n_points == 0
            assert sp1.points_dist_type == "g"
    #
    sp_l1_list = ssp_r.shuffle_sp_list(
        q.SelectedPointsComplexD, sp2_list, is_reverse=True
    )
    sig_l1_list = get_f_list_sig(sp_l1_list, rs, 3)
    q.json_results_append("sig sp_l1_list", sig_l1_list, 1e-12)
    assert np.all(sig_l1_list == sig_l_list)

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
    q.json_results_append("get_data_sig_arr(sp_l,rs,2)", sig_l, 1e-12)
    sp_s = q.SelectedPointsComplexD(sp_l, ssp)
    assert sp_s.psel is ssp.psel_dst_list[0]
    sig_s = q.get_data_sig_arr(sp_s, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp_s,rs,2)", sig_s, 1e-12)
    psel_l1 = q.PointsSelection(psel_s, ssp, True)
    assert psel_l1 == psel_l
    sp_l1 = q.SelectedPointsComplexD(sp_s, ssp, True)
    assert sp_l1.psel is ssp.psel_src_list[0]
    sig_l1 = q.get_data_sig_arr(sp_l1, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp_l1,rs,2)", sig_l1, 1e-12)
    assert np.all(sig_l1 == sig_l)
    #
    sp_lc = q.SelectedPointsChar()
    sp_l1.swap_cast(sp_lc)
    sp_sc = q.SelectedPointsChar(sp_lc, ssp)
    sp_l1.swap_cast(sp_sc)
    sp_s1 = sp_l1
    sig_s1 = q.get_data_sig_arr(sp_s1, rs.split("sig"), 2)
    q.json_results_append("get_data_sig_arr(sp_s1,rs,2)", sig_s1, 1e-12)
    assert np.all(sig_s1 == sig_s)
    #
    psel_l_list = [psel_l.copy() for i in range(3)]
    geo_list = [geo.copy() for i in range(3)]
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
    q.json_results_append("sig sp_l_list", get_f_list_sig(sp_l_list, rs, 3), 1e-12)
    #
    sp_s_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_l_list)
    q.json_results_append("sig sp_s_list", get_f_list_sig(sp_s_list, rs, 3), 1e-12)
    sp_ss_list = ssp.shuffle_sp_list(
        q.SelectedPointsComplexD, sp_s_list, is_reverse=True
    )
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
    psel_l_list = [psel_l.copy() for i in range(3)]
    geo_list = [geo.copy() for i in range(3)]
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
    q.json_results_append("sig sp_l_list", get_f_list_sig(sp_l_list, rs, 3), 1e-12)
    #
    sp_s_list = ssp.shuffle_sp_list(q.SelectedPointsComplexD, sp_l_list)
    q.json_results_append("sig sp_s_list", get_f_list_sig(sp_s_list, rs, 3), 1e-12)
    sp_ss_list = ssp.shuffle_sp_list(
        q.SelectedPointsComplexD, sp_s_list, is_reverse=True
    )
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
    q.json_results_append("get_data_sig_arr(sp_l,rs,2)", sig_l, 1e-12)
    #
    sp_s_list = ssp.shuffle_sp_list(
        q.SelectedPointsComplexD,
        [
            sp_l,
        ],
    )
    q.json_results_append("sig sp_s_list", get_f_list_sig(sp_s_list, rs, 3), 1e-12)
    sp_ss_list = ssp.shuffle_sp_list(
        q.SelectedPointsComplexD, sp_s_list, is_reverse=True
    )
    assert np.all(
        get_f_list_sig(sp_ss_list, rs, 3)
        == get_f_list_sig(
            [
                sp_l,
            ],
            rs,
            3,
        )
    )

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
    t_slice_size_node = q.Coordinate(
        [
            1,
            1,
            1,
            total_site[3],
        ]
    )
    #
    num_field = 3
    #
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split("gf"), 0.5, 2)
    gf_sig = q.get_data_sig_arr(gf, rs, 3)
    q.json_results_append("gf sig", gf_sig, 1e-10)
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
    q.json_results_append("sig vec_list", vec_sig, 1e-12)
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
    q.json_results_append(
        f"hash_sha256(ssp1.psel_dst_list)={q.hash_sha256(ssp1.psel_dst_list)}"
    )
    t_list = [psel[0][3].item() for psel in ssp1.psel_dst_list]
    s_geo_list = [
        q.Geometry(
            q.Coordinate(
                [
                    0,
                    0,
                    0,
                    t,
                ]
            ),
            t_slice_size_node,
            t_slice_node_site,
        )
        for t in t_list
    ]
    all_t_list_list = q.get_comm().allgather(t_list)
    q.json_results_append(f"all_t_list_list={all_t_list_list}")
    #
    gf.swap_sp_cast(spc, geo)
    s_spc_list = ssp1.shuffle_list(
        [
            spc,
        ]
    )
    s_gf_list = []
    for s_spc, s_geo in zip(s_spc_list, s_geo_list):
        s_gf = q.GaugeField()
        s_gf.swap_sp_cast(s_spc, s_geo)
        s_gf_list.append(s_gf)
    q.json_results_append("sig s_gf_list", get_f_list_sig(s_gf_list, rs, 3), 1e-10)
    for s_gf, s_spc, s_geo in zip(s_gf_list, s_spc_list, s_geo_list):
        s_gf.swap_sp_cast(s_spc, s_geo)
    [
        spc,
    ] = ssp1.shuffle_list(s_spc_list, is_reverse=True)
    gf.swap_sp_cast(spc, geo)
    assert np.all(gf_sig == q.get_data_sig_arr(gf, rs, 3))
    #
    ssp2 = q.SelectedShufflePlan("t_slice_from_l", psel_list, geo_list)
    assert ssp2.psel_src_list is psel_list
    q.json_results_append(
        f"hash_sha256(ssp2.psel_dst_list)={q.hash_sha256(ssp2.psel_dst_list)}"
    )
    vt_list = [psel[0][3].item() for psel in ssp2.psel_dst_list]
    vs_geo_list = [
        q.Geometry(
            q.Coordinate(
                [
                    0,
                    0,
                    0,
                    t,
                ]
            ),
            t_slice_size_node,
            t_slice_node_site,
        )
        for t in vt_list
    ]
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
    q.json_results_append("sig s_vec_list", get_f_list_sig(s_vec_list, rs, 3), 1e-10)
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
    gf.set_rand(rs.split("gf"), 0.5, 2)
    gf_sig = q.get_data_sig_arr(gf, rs, 3)
    q.json_results_append("gf sig", gf_sig, 1e-10)
    #
    ff_list = []
    for id_field in range(num_field):
        ff = q.FermionField4d(geo)
        ff.set_rand(rs.split(f"ff_list[{id_field}]"))
        ff_list.append(ff)
    ff_list_sig = get_f_list_sig(ff_list, rs, 3)
    q.json_results_append("sig ff_list", ff_list_sig, 1e-12)
    #
    s_ff_list = q.prop_spatial_smear(ff_list, gf, coef, step)
    s_ff_list_sig = get_f_list_sig(s_ff_list, rs, 3)
    q.json_results_append("sig s_ff_list", s_ff_list_sig, 1e-12)
    #
    mom = q.CoordinateD(
        [
            0.1,
            -0.2,
            0.3,
            0.5,
        ]
    )
    #
    ms_ff_list = q.prop_spatial_smear(ff_list, gf, coef, step, mom, chunk_size=2)
    ms_ff_list_sig = get_f_list_sig(ms_ff_list, rs, 3)
    q.json_results_append("sig ms_ff_list", ms_ff_list_sig, 1e-12)
    #
    if multiplicity == 1:
        prop = q.Prop(geo)
        prop.set_rand(rs.split("prop"))
        prop_sig = q.get_data_sig_arr(prop, rs, 3)
        q.json_results_append("prop sig", prop_sig, 1e-10)
        ss_prop = q.prop_spatial_smear(prop, gf, coef, step, mom)
        ss_prop_sig = q.get_data_sig_arr(ss_prop, rs, 3)
        q.json_results_append("ss_prop sig", ss_prop_sig, 1e-10)
        ss_prop = q.prop_smear(prop, gf, coef, step, mom)
        ss_prop_sig = q.get_data_sig_arr(ss_prop, rs, 3)
        q.json_results_append("ss_prop sig", ss_prop_sig, 1e-10)
        ss_prop = q.prop_smear(prop, gf, coef, step, mom, mode_smear=0)
        ss_prop_sig = q.get_data_sig_arr(ss_prop, rs, 3)
        q.json_results_append("ss_prop sig", ss_prop_sig, 1e-10)

for total_site in total_site_list:
    for multiplicity in multiplicity_list:
        for seed in range(1):
            selected_shuffle_l_from_g(total_site, multiplicity, seed)
            selected_shuffle_r_from_l(total_site, multiplicity, seed)
            selected_shuffle_t_slice_from_l(total_site, multiplicity, seed)
            selected_shuffle_t_slice_from_f(total_site, multiplicity, seed)
            test_prop_spatial_smear(total_site, multiplicity, seed)

q.timer_display()
q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
