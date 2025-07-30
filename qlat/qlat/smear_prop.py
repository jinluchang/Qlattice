__all__ = [
        'prop_spatial_smear',
        ]

import numpy as np

class q:
    from qlat_utils import (
            timer,
            get_chunk_list,
            CoordinateD,
            )
    from .c import (
            GaugeField,
            Prop,
            FermionField4d,
            PointsSelection,
            SelectedShufflePlan,
            prop_spatial_smear_no_comm,
            mk_ff_list_from_prop,
            mk_prop_from_ff_list,
            )

@q.timer
def prop_spatial_smear(ff_list, gf, coef, step, mom=None, *, chunk_size=12):
    """
    Perform spatial smear for `ff_list`, a list of `FermionField4d`.
    Return new `ff_list` after smearing.
    Original `ff_list` should not be modified.
    #
    `ff_list` can also be a `Prop`.
    In this case, the function will also return a `Prop`.
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
    #
    if isinstance(ff_list, q.Prop):
        prop = ff_list
        ff_list = q.mk_ff_list_from_prop(prop)
        ss_ff_list = prop_spatial_smear(ff_list, gf, coef, step, mom, chunk_size=chunk_size)
        ss_prop = q.mk_prop_from_ff_list(ss_ff_list)
        return ss_prop
    #
    assert isinstance(gf, q.GaugeField)
    assert isinstance(coef, float)
    assert isinstance(step, int)
    assert step >= 0
    if mom is None:
        mom = q.CoordinateD()
    assert isinstance(mom, q.CoordinateD)
    #
    assert isinstance(ff_list, list)
    for ff in ff_list:
        assert isinstance(ff, q.FermionField4d)
    #
    geo = gf.geo
    for ff in ff_list:
        assert geo == ff.geo
    #
    if step == 0:
        return
    #
    chunk_ff_list_list = q.get_chunk_list(ff_list, chunk_size=chunk_size)
    #
    ss_ff_list = []
    #
    plan = dict()
    #
    for chunk_ff_list in chunk_ff_list_list:
        num_field = len(chunk_ff_list)
        if plan.get("num_field") != num_field:
            plan = prop_spatial_smear_chunk_planner(gf, num_field)
        assert plan["num_field"] == len(chunk_ff_list)
        ss_chunk_ff_list = prop_spatial_smear_chunk(chunk_ff_list, plan, coef, step, mom)
        ss_ff_list += ss_chunk_ff_list
    #
    return ss_ff_list

@q.timer
def prop_spatial_smear_chunk_planner(gf, num_field):
    """
    return plan
    #
    plan = dict(s_gf_list=s_gf_list, ssp=ssp)
    """
    geo = gf.geo
    total_site = geo.total_site
    psel = q.PointsSelection(geo)
    psel_list = [ psel.copy() for i in range(num_field) ]
    geo_list = [ geo.copy() for i in range(num_field) ]
    #
    ssp_gf = q.SelectedShufflePlan("dist_t_slice_from_l", psel, geo, num_field)
    ssp = q.SelectedShufflePlan("t_slice_from_l", psel_list, geo_list)
    #
    s_gf_list = ssp_gf.shuffle_sp_list(q.GaugeField, [ gf, ])
    #
    plan = dict(
            num_field=num_field,
            s_gf_list=s_gf_list,
            ssp=ssp,
            )
    return plan

@q.timer
def prop_spatial_smear_chunk(ff_list, plan, coef, step, mom):
    """
    return ss_ff_list
    which is the smeared FermionField4d fields.
    """
    num_field = plan['num_field']
    s_gf_list = plan['s_gf_list']
    ssp = plan['ssp']
    #
    assert len(ff_list) == num_field
    #
    s_ff_list = ssp.shuffle_sp_list(q.FermionField4d, ff_list)
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
    ss_ff_list = ssp.shuffle_sp_list(q.FermionField4d, s_ff_list, is_reverse=True)
    return ss_ff_list
