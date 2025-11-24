#!/usr/bin/env python3

import sys
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

rs = q.RngState("seed")

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.5, 10)
q.json_results_append(f"gf-init plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf-init", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)
gf0 = gf.copy()

gf @= gf0
q.gf_block_stout_smear(gf, q.Coordinate(), 0.1)
q.json_results_append(f"gf_block_stout_smear after plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
q.gf_block_stout_smear(gf, q.Coordinate([ 2, 2, 2, 2, ]), 0.1)
q.json_results_append(f"gf_block_stout_smear after block plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after block", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
q.gf_block_stout_smear(gf, q.Coordinate([ 4, 4, 4, 4, ]), 0.1)
q.json_results_append(f"gf_block_stout_smear after block 4 plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after block 4", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
gf = q.field_shift(gf, q.Coordinate([ 1, 2, 3, 4, ]))
q.json_results_append(f"field_shift plaq", gf.plaq(), 1e-8)
q.json_results_append(f"field_shift", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
q.gf_local_stout_smear(gf, q.Coordinate([ 4, 4, 4, 4, ]), 0.1)
q.json_results_append(f"gf_local_stout_smear after block 4 plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_local_stout_smear after block 4", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
q.gf_local_stout_smear(gf, q.Coordinate([ 2, 2, 2, 2, ]), 0.1)
q.json_results_append(f"gf_local_stout_smear after block 2 plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_local_stout_smear after block 2", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

fc = q.FieldChar(geo, 1)
fc[:] = np.arange(len(fc[:]))[:, None]
fc_init = fc[:].copy()
fc_list = q.shuffle_field_char(fc, q.Coordinate([ 2, 2, 2, 4, ]))
q.json_results_append(f"shuffle_field_char {len(fc_list)}")
fc[:] = 0
q.shuffle_field_char_back(fc, fc_list, q.Coordinate([ 2, 2, 2, 4, ]))
assert np.all(fc[:] == fc_init)

gf @= gf0
gf_list = q.shuffle_field(gf, q.Coordinate([ 2, 2, 2, 4, ]))
q.json_results_append(f"shuffle_field {len(gf_list)}")
for gf_local in gf_list:
    q.gf_local_stout_smear(gf_local, q.Coordinate([ 2, 2, 2, 2, ]), 0.1)
q.shuffle_field_back(gf, gf_list, q.Coordinate([ 2, 2, 2, 4, ]))
q.json_results_append(f"gf_local_stout_smear plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_local_stout_smear", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gt = q.GaugeTransform(geo)
gt.set_rand(rs.split("gt-init"), 0.5, 2)
q.json_results_append(f"gt-init", q.get_data_sig_arr(gt, q.RngState(), 3), 1e-8)
gt_inv = gt.inv()
q.json_results_append(f"gt_inv", q.get_data_sig_arr(gt_inv, q.RngState(), 3), 1e-8)

gf_gt = gt * gf
q.json_results_append(f"gf_gt", q.get_data_sig_arr(gf_gt, q.RngState(), 3), 1e-8)

gt1 = gt * gt_inv
q.json_results_append(f"gt1", q.get_data_sig_arr(gt1, q.RngState(), 3), 1e-8)

prop = q.Prop(geo)
prop.set_rand_g(rs.split("prop-init"), 0.0, 0.5)
q.json_results_append(f"prop", q.get_data_sig_arr(prop, q.RngState(), 3), 1e-8)

prop_gt = gt * prop
q.json_results_append(f"prop_gt", q.get_data_sig_arr(prop_gt, q.RngState(), 3), 1e-8)

psel = q.PointsSelection()
psel.set_rand(geo.total_site, 16, rs.split("psel"))
sp_prop = q.PselProp(psel)
sp_prop @= prop
q.json_results_append(f"sp_prop", q.get_data_sig_arr(sp_prop, q.RngState(), 3), 1e-8)

sp_prop_gt = gt * sp_prop
q.json_results_append(f"sp_prop_gt", q.get_data_sig_arr(sp_prop_gt, q.RngState(), 3), 1e-8)
sp_prop_gt.set_zero()
sp_prop_gt @= prop_gt
q.json_results_append(f"sp_prop_gt 2", q.get_data_sig_arr(sp_prop_gt, q.RngState(), 3), 1e-8)

fsel = q.FieldSelection(psel)
s_prop = q.SelProp(fsel)
s_prop @= prop
q.json_results_append(f"s_prop", q.get_data_sig_arr(s_prop, q.RngState(), 3), 1e-8)

s_prop_gt = gt * s_prop
q.json_results_append(f"s_prop_gt", q.get_data_sig_arr(s_prop_gt, q.RngState(), 3), 1e-8)

s_prop_gt.set_zero()
s_prop_gt @= prop_gt
q.json_results_append(f"s_prop_gt 2", q.get_data_sig_arr(s_prop_gt, q.RngState(), 3), 1e-8)

ff4d = q.FermionField4d(geo)
ff4d.set_rand_g(rs.split("ff-init"), 0.0, 1.0)
q.json_results_append(f"ff4d", q.get_data_sig_arr(ff4d, q.RngState(), 3), 1e-8)

ff4d_gt = gt * ff4d
q.json_results_append(f"ff4d_gt", q.get_data_sig_arr(ff4d_gt, q.RngState(), 3), 1e-8)

gf @= gf0
gf_list = q.shuffle_field(gf, q.Coordinate([ 2, 2, 2, 4, ]))
geo_list = [ gf_local.geo for gf_local in gf_list ]
rs_gt_tree = rs.split("gt_tree")
f_dir_list = [ q.mk_local_tree_gauge_f_dir(geo, q.Coordinate([ 2, 2, 2, 2, ]), rs_gt_tree) for geo in geo_list ]
gt_inv_list = []
q.json_results_append(f"shuffle_field {len(gf_list)}")
for gf_local, f_dir in zip(gf_list, f_dir_list):
    for step in range(1):
        q.gf_local_stout_smear(gf_local, q.Coordinate([ 2, 2, 2, 2, ]), 0.1)
    gt_inv = q.gt_local_tree_gauge(gf_local, f_dir, 4)
    gt_inv_list.append(gt_inv)
q.shuffle_field_back(gf, gf_list, q.Coordinate([ 2, 2, 2, 4, ]))
q.json_results_append(f"gf_local_stout_smear plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_local_stout_smear", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)
q.json_results_append(f"gf_local_stout_smear trace", gf.link_trace(), 1e-8)

gt_inv = q.GaugeTransform(geo)
q.shuffle_field_back(gt_inv, gt_inv_list, q.Coordinate([ 2, 2, 2, 4, ]))
q.json_results_append(f"gt_inv", q.get_data_sig_arr(gt_inv, q.RngState(), 3), 1e-8)

gt = gt_inv.inv()
gf @= gt * gf
q.json_results_append(f"gf_local_stout_smear gt plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_local_stout_smear gt", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)
q.json_results_append(f"gf_local_stout_smear gt trace", gf.link_trace(), 1e-8)

gf @= gf0
f_dir_list = None
gt_inv, f_dir_list = q.gt_block_tree_gauge(
    gf,
    block_site=q.Coordinate([ 2, 2, 2, 2, ]),
    stout_smear_step_size=0.1,
    num_smear_step=1,
    f_dir_list=f_dir_list,
    rs_f_dir=rs.split("gt_tree"),
    )
q.json_results_append(f"gt_block_tree_gauge gt_inv", q.get_data_sig_arr(gt_inv, q.RngState(), 3), 1e-8)
q.json_results_append(f"gf", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
f_dir_list = None
gt_inv, f_dir_list = q.gt_block_tree_gauge(
    gf,
    block_site=q.Coordinate([ 2, 2, 2, 2, ]),
    new_size_node=q.Coordinate([ 1, 1, 1, 2, ]),
    stout_smear_step_size=0.1,
    num_smear_step=1,
    f_dir_list=f_dir_list,
    rs_f_dir=rs.split("gt_tree"),
    )
q.json_results_append(f"gt_block_tree_gauge gt_inv", q.get_data_sig_arr(gt_inv, q.RngState(), 3), 1e-8)
q.json_results_append(f"gf", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gt = gt_inv.inv()
gf @= gt * gf
q.json_results_append(f"gf_local_stout_smear gt_block_tree_gauge plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_local_stout_smear gt_block_tree_gauge", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)
q.json_results_append(f"gf_local_stout_smear gt_block_tree_gauge trace", gf.link_trace(), 1e-8)

gf @= gf0

f_dir_list = None
gt_inv, f_dir_list = q.gt_block_tree_gauge(
    gf,
    block_site=q.Coordinate([ 4, 4, 4, 4, ]),
    f_dir_list=f_dir_list,
    rs_f_dir=rs.split("gt_tree"),
    )
q.json_results_append(f"gt_block_tree_gauge 4 gt_inv", q.get_data_sig_arr(gt_inv, q.RngState(), 3), 1e-8)
q.json_results_append(f"gf after gt_block_tree_gauge", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf = gt_inv.inv() * gf
q.json_results_append(f"gf gt transformed", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gt_inv, f_dir_list = q.gt_block_tree_gauge(
    gf,
    block_site=q.Coordinate([ 4, 4, 4, 4, ]),
    f_dir_list=f_dir_list,
    rs_f_dir=rs.split("gt_tree"),
    )
q.json_results_append(f"gt_block_tree_gauge 4 gt_inv 2nd", q.get_data_sig_arr(gt_inv, q.RngState(), 3), 1e-8)
q.json_results_append(f"gf after gt_block_tree_gauge 2nd", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf = gt_inv.inv() * gf
q.json_results_append(f"gf gt transformed 2nd", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gm = q.GaugeMomentum(geo)
gm.set_rand(rs.split("gm-init"), 1.0)
q.json_results_append(f"gm", q.get_data_sig_arr(gm, q.RngState(), 3), 1e-8)

basis = q.FieldRealD()
q.set_basis_from_anti_hermitian_matrix(basis, gm)
q.json_results_append(f"basis", q.get_data_sig_arr(basis, q.RngState(), 3), 1e-8)

gm.set_unit()
q.set_anti_hermitian_matrix_from_basis(gm, basis)
q.json_results_append(f"gm", q.get_data_sig_arr(gm, q.RngState(), 3), 1e-8)

q.timer_display()
q.check_log_json(__file__, check_eps=1e-5)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
