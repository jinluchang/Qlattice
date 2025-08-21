r"""

# Plaq dependent flow

Standard Wilson flow action:
$$
S_\mathrm{Wilson} = \frac{\beta}{2}\sum_{x,\mu,\nu} \Big(1 - \frac{1}{3}\mathrm{Re}\mathrm{Tr} U_{\mu,\nu}\Big)
$$
with $\beta = 3$.

Modified flow to make the slope depends on the plaq value:
$$
S_f = -\frac{\beta}{2}\sum_{x,\mu,\nu} f\big(\frac{1}{3}\mathrm{Re}\mathrm{Tr} U_{\mu,\nu}\big)
$$

We can choose the function $f$ that smooths gauge field and therefore fix the total topological charge:
$$
\frac{d}{dp}f_\mathrm{Freeze}(p)
=
1 - p
$$
Note that $f_\mathrm{Freeze}$ suppresses plaq significantly deviate from $1$ and therefore prevent topological charge tunnelling during flow.

We can also choose the function $f$ that shrinks instanton:
$$
\frac{d}{dp}f_\mathrm{Shrink}(p)
=
\frac{\epsilon}{1 - p + \epsilon} + b
$$
with a possible choice of parameters be $\epsilon = 0.005$ and $b=0.5$.
Note that $f_\mathrm{Shrink}$ tends to shrink the size of instanton therefore enhance topological charge tunnelling during flow.

"""

### --------------------------------------

import numpy as np

class q:
    from qlat_utils import (
            timer,
            Coordinate,
            smod_coordinate_d,
            get_fname,
            RngState,
            parallel_map,
            random_permute,
            displayln_info,
            get_num_node,
            is_test,
            json_results_append,
            )
    from .c import (
            Geometry,
            PointsSelection,
            FieldComplexD,
            FieldRealD,
            gf_plaq_flow_force,
            gf_evolve,
            gf_plaq_field,
            gf_topology_field,
            gf_topology,
            Coordinate,
            CoordinateD,
            PointsSelection,
            gf_reduce_half,
            )
    from .field_analysis import (
            smear_field,
            sphere_sum_field,
            )
    from .mpi_utils import (
            get_mpi_chunk,
            get_comm,
            )

### --------------------------------------

@q.timer
def gf_flow_topo(gf, step_size, tag=None):
    """
    Modify `gf` in place.
    Default Wilson flow with Euler integrator.
    tag in [ None, "Shrink", "Freeze", ]
    """
    geo = gf.geo
    plaq_factor = q.FieldRealD(geo, 6)
    if tag is None:
        plaq_factor.set_unit()
    else:
        f_plaq = q.gf_plaq_field(gf)
        plaq_min = np.min(f_plaq.glb_min()[:])
        plaq_max = np.max(f_plaq.glb_max()[:])
        if tag == "Shrink":
            eps = 0.005
            base = 0.5
            norm = eps / (1 - plaq_max + eps) + base
            plaq_factor[:] = eps / (1 - f_plaq[:] + eps) + base
            plaq_factor[:] *= 1 / norm
        elif tag == "Freeze":
            norm = 1 - plaq_min
            plaq_factor[:] = 1 - f_plaq[:]
            plaq_factor[:] *= 1 / norm
        else:
            fname = q.get_fname()
            raise Exception(f"{fname}: tag={tag}")
    gm_force = q.gf_plaq_flow_force(gf, plaq_factor)
    q.gf_evolve(gf, gm_force, step_size)

@q.timer
def mk_plaq_xg_arr(geo):
    n_plaq = 6
    xg_arr = geo.xg_arr()
    (n_points, n_dim,) = xg_arr.shape
    assert n_dim == 4
    assert n_points == geo.local_volume
    plaq_xg_arr = np.zeros((n_points, n_plaq, n_dim,), dtype=np.float64)
    plaq_xg_arr[:, :, :] = xg_arr[:, None, :]
    plaq_xg_arr[:, 0, :] += [ 0.5, 0.5, 0, 0, ]
    plaq_xg_arr[:, 1, :] += [ 0.5, 0, 0.5, 0, ]
    plaq_xg_arr[:, 2, :] += [ 0.5, 0, 0, 0.5, ]
    plaq_xg_arr[:, 3, :] += [ 0, 0.5, 0.5, 0, ]
    plaq_xg_arr[:, 4, :] += [ 0, 0.5, 0, 0.5, ]
    plaq_xg_arr[:, 5, :] += [ 0, 0, 0.5, 0.5, ]
    # print(plaq_xg_arr.shape)
    return plaq_xg_arr

@q.timer
def get_extreme_plaq_xg_list(plaq_xg_arr, f_plaq, threshold):
    """
    return extreme_plaq_xg_list
    extreme_plaq_xg_list = [ (plaq, xg,), ... ]
    #
    e.g.
    plaq_xg_arr = mk_plaq_xg_arr(geo)
    f_plaq = q.gf_plaq_field(gf)
    threshold = 0.1
    extreme_plaq_xg_list = get_extreme_plaq_xg_list(plaq_xg_arr, f_plaq, threshold)
    """
    sel = f_plaq[:] < (1 - threshold)
    plaq_arr = f_plaq[sel]
    extreme_xg_arr = plaq_xg_arr[sel]
    extreme_plaq_xg_list = list(zip(plaq_arr.tolist(), extreme_xg_arr.tolist()))
    all_extreme_plaq_xg_list = q.get_comm().allgather(extreme_plaq_xg_list)
    extreme_plaq_xg_list = []
    for sub in all_extreme_plaq_xg_list:
        extreme_plaq_xg_list += sub
    extreme_plaq_xg_list.sort()
    return extreme_plaq_xg_list

def point_d_dis_sqr(x, y, total_site):
    return q.smod_coordinate_d(x - y, total_site).sqr()

@q.timer
def get_group_extreme_plaq_xg_list(extreme_plaq_xg_list, total_site, dis_sqr_limit):
    """
    Group points with distance square less than `dis_sqr_limit`.
    #
    return group_extreme_plaq_xg_list
    group_extreme_plaq_xg_list = [
        [ (plaq, xg, dis_sq,), ... ],
        ...
    ]
    #
    e.g.
    extreme_plaq_xg_list = get_extreme_plaq_xg_list(plaq_xg_arr, f_plaq, threshold)
    total_site = q.Coordinate([ 4, 4, 4, 8, ])
    dis_sqr_limit = 5
    group_extreme_plaq_xg_list = get_group_extreme_plaq_xg_list(extreme_plaq_xg_list, total_site, dis_sqr_limit)
    """
    total_site = q.CoordinateD(total_site)
    group_extreme_plaq_xg_list = []
    for plaq, xg in extreme_plaq_xg_list:
        x1 = q.CoordinateD(xg)
        b = False
        for l in group_extreme_plaq_xg_list:
            assert len(l) >= 1
            x2 = q.CoordinateD(l[0][1])
            dis_sq = point_d_dis_sqr(x1, x2, total_site)
            if dis_sq < dis_sqr_limit:
                l.append((plaq, xg, dis_sq))
                b = True
                break
        if not b:
            l = [ (plaq, xg, 0,), ]
            group_extreme_plaq_xg_list.append(l)
    return group_extreme_plaq_xg_list

@q.timer
def process_inst_list(inst_list):
    """
    inst_list[inst_idx] = [
        dict(flow_time, current_spacing, total_site, o_xg_d, o_xg, xg_d, xg, plaq, s_topo),
        ...
    ]
    #
    p_inst_list = process_inst_list(inst_list)
    #
    p_inst_list[inst_idx] = dict(
        # original `inst_idx` for `inst_list`
        inst_idx,
        inst_raw = inst_list[inst_idx],
        #
        # Values are usually taken at minimum `plaq`
        o_xg_d,
        o_total_site,
        xg_d,
        total_site,
        current_spacing,
        plaq,
        flow_time,
        num_plaq,
        dis_sqr_max,
        #
        flow_time_start,
        flow_time_end,
        dis_sqr_max_max,
        c_dis_sqr_max,
        #
        plaq_list,
        s_topo_list,
        num_plaq_list,
        dis_sqr_avg_list,
        dis_sqr_max_list,
        c_dis_sqr_list,
        delta_s_topo,
    )
    """
    p_inst_list = []
    for inst_idx, inst in enumerate(inst_list):
        s_topo_list = []
        plaq_list = []
        num_plaq_list = []
        dis_sqr_avg_list = []
        dis_sqr_max_list = []
        c_dis_sqr_list = []
        plaq = 1.0
        for info_dict in inst:
            plaq_list.append(info_dict["plaq"])
            num_plaq_list.append(info_dict["num_plaq"])
            dis_sqr_avg_list.append(info_dict["dis_sqr_avg"])
            dis_sqr_max_list.append(info_dict["dis_sqr_max"])
            if info_dict["s_topo"] is not None:
                s_topo_list.append(info_dict["s_topo"])
            if info_dict["plaq"] < plaq:
                plaq = info_dict["plaq"]
                flow_time = info_dict["flow_time"]
                s_topo = info_dict["s_topo"]
                num_plaq = info_dict["num_plaq"]
                dis_sqr_avg = info_dict["dis_sqr_avg"]
                dis_sqr_max = info_dict["dis_sqr_max"]
                o_xg_d = info_dict["o_xg_d"]
                xg_d = info_dict["xg_d"]
        for info_dict in inst:
            total_site_d = q.CoordinateD(info_dict["total_site"])
            c_dis_sqr = point_d_dis_sqr(q.CoordinateD(xg_d), q.CoordinateD(info_dict["xg_d"]), total_site_d)
            c_dis_sqr_list.append(c_dis_sqr)
        s_topo_arr = np.array(s_topo_list, dtype=np.float64)
        s_topo_list = s_topo_arr.T.tolist()
        delta_s_topo = None
        if len(s_topo_arr) > 0:
            delta_s_topo = (s_topo_arr[0] - s_topo_arr[-1]).tolist()
        total_site = inst[0]["total_site"]
        current_spacing = inst[0]["current_spacing"]
        topo_sphere_sum_radius_list = inst[0]["topo_sphere_sum_radius_list"]
        dis_sqr_max_max = np.max(dis_sqr_max_list).item()
        c_dis_sqr_max = np.max(c_dis_sqr_list).item()
        p_inst = dict()
        p_inst["inst_idx"] = inst_idx
        p_inst["inst_raw"] = inst
        p_inst["o_xg_d"] = o_xg_d
        p_inst["o_total_site"] = [ current_spacing * c for c in total_site ]
        p_inst["xg_d"] = xg_d
        p_inst["total_site"] = total_site
        p_inst["current_spacing"] = current_spacing
        p_inst["plaq"] = plaq
        p_inst["flow_time"] = flow_time
        p_inst["flow_time_start"] = inst[0]["flow_time"]
        p_inst["flow_time_end"] = inst[-1]["flow_time"]
        p_inst["num_plaq"] = num_plaq
        p_inst["dis_sqr_avg"] = dis_sqr_avg
        p_inst["dis_sqr_max"] = dis_sqr_max
        p_inst["dis_sqr_max_max"] = dis_sqr_max_max
        p_inst["c_dis_sqr_max"] = c_dis_sqr_max
        p_inst["delta_s_topo"] = delta_s_topo
        p_inst["plaq_list"] = plaq_list
        p_inst["s_topo_list"] = s_topo_list
        p_inst["topo_sphere_sum_radius_list"] = topo_sphere_sum_radius_list
        p_inst["num_plaq_list"] = num_plaq_list
        p_inst["dis_sqr_avg_list"] = dis_sqr_avg_list
        p_inst["dis_sqr_max_list"] = dis_sqr_max_list
        p_inst["c_dis_sqr_list"] = c_dis_sqr_list
        p_inst_list.append(p_inst)
    def mk_closest_inst_info(inst_dis_sqr, p_inst):
        """
        e.g.
        mk_closest_inst_info(inst_dis_sqr, p_inst_2)
        """
        closest_inst_info = dict()
        closest_inst_info["inst_dis_sqr"] = inst_dis_sqr
        closest_inst_info["current_spacing"] = p_inst["current_spacing"]
        closest_inst_info["inst_idx"] = p_inst["inst_idx"]
        closest_inst_info["o_xg_d"] = p_inst["o_xg_d"]
        closest_inst_info["plaq"] = p_inst["plaq"]
        closest_inst_info["flow_time"] = p_inst["flow_time"]
        closest_inst_info["flow_time_start"] = p_inst["flow_time_start"]
        closest_inst_info["flow_time_end"] = p_inst["flow_time_end"]
        closest_inst_info["dis_sqr_max_max"] = p_inst["dis_sqr_max_max"]
        closest_inst_info["c_dis_sqr_max"] = p_inst["c_dis_sqr_max"]
        closest_inst_info["delta_s_topo"] = p_inst["delta_s_topo"]
        return closest_inst_info
    for p_inst in p_inst_list:
        inst_idx = p_inst["inst_idx"]
        o_xg_d = p_inst["o_xg_d"]
        o_xg_d = q.CoordinateD(o_xg_d)
        o_total_site = p_inst["o_total_site"]
        o_total_site_d = q.CoordinateD(o_total_site)
        current_spacing = p_inst["current_spacing"]
        flow_time = p_inst["flow_time"]
        flow_time_start = p_inst["flow_time_start"]
        flow_time_end = p_inst["flow_time_end"]
        closest_inst_dis_sqr = None
        closest_inst_info = None
        closest_sim_inst_dis_sqr = None
        closest_sim_inst_info = None
        for p_inst_2 in p_inst_list:
            inst_idx_2 = p_inst_2["inst_idx"]
            if inst_idx_2 == inst_idx:
                continue
            o_xg_d_2 = p_inst_2["o_xg_d"]
            o_xg_d_2 = q.CoordinateD(o_xg_d_2)
            o_total_site_2 = p_inst_2["o_total_site"]
            assert o_total_site_2 == o_total_site
            current_spacing_2 = p_inst_2["current_spacing"]
            flow_time_start_2 = p_inst_2["flow_time_start"]
            flow_time_end_2 = p_inst_2["flow_time_end"]
            inst_dis_sqr = point_d_dis_sqr(o_xg_d, o_xg_d_2, o_total_site_d)
            if (closest_inst_dis_sqr is None) or (inst_dis_sqr < closest_inst_dis_sqr):
                closest_inst_dis_sqr = inst_dis_sqr
                closest_inst_info = mk_closest_inst_info(inst_dis_sqr, p_inst_2)
            if ((flow_time_start <= flow_time_start_2 <= flow_time_end)
                or (flow_time_start <= flow_time_end_2 <= flow_time_end)
                or (flow_time_start_2 <= flow_time_start <= flow_time_end_2)
                or (flow_time_start_2 <= flow_time_end <= flow_time_end_2)
                ):
                if (closest_sim_inst_dis_sqr is None) or (inst_dis_sqr < closest_sim_inst_dis_sqr):
                    assert current_spacing_2 == current_spacing
                    closest_sim_inst_dis_sqr = inst_dis_sqr
                    closest_sim_inst_info = mk_closest_inst_info(inst_dis_sqr, p_inst_2)
        p_inst["closest_sim_inst_info"] = closest_sim_inst_info
        p_inst["closest_inst_info"] = closest_inst_info
    for p_inst in p_inst_list:
        delta_s_topo = p_inst["delta_s_topo"]
        estimate_topo_charge = 0
        p_inst["estimate_topo_charge"] = estimate_topo_charge
        if delta_s_topo is None:
            continue
        delta_s_topo = np.array(delta_s_topo, dtype=np.float64)
        if len(delta_s_topo) > 1:
            estimate_topo_charge = (delta_s_topo[0] - (delta_s_topo[1] - 1.5 * delta_s_topo[0]) / 3) / 0.44
        elif len(delta_s_topo) > 0:
            estimate_topo_charge = delta_s_topo[0] / 0.44
        p_inst["estimate_topo_charge"] = estimate_topo_charge
        #
        # # If there is another instanton very close, consider the two instantons together. (Disable now)
        # p_inst["estimate_topo_charge_orig"] = None
        # inst_idx = p_inst["inst_idx"]
        # current_spacing = p_inst["current_spacing"]
        # closest_sim_inst_info = p_inst["closest_sim_inst_info"]
        # if closest_sim_inst_info is None:
        #     continue
        # inst_dis_sqr_2 = closest_sim_inst_info["inst_dis_sqr"]
        # pairwise_inst_dis_sqr_max_limit = 4**2
        # pairwise_inst_dis_sqr_min_limit = 5
        # if inst_dis_sqr_2 >= pairwise_inst_dis_sqr_max_limit * current_spacing**2:
        #     continue
        # if inst_dis_sqr_2 <= pairwise_inst_dis_sqr_min_limit * current_spacing**2:
        #     continue
        # inst_idx_2 = closest_sim_inst_info["inst_idx"]
        # closest_sim_inst_info_2 = p_inst_list[inst_idx_2]["closest_sim_inst_info"]
        # assert closest_sim_inst_info_2 is not None
        # if closest_sim_inst_info_2["inst_idx"] != inst_idx:
        #     continue
        # delta_s_topo_2 = closest_sim_inst_info["delta_s_topo"]
        # if delta_s_topo_2 is None:
        #     continue
        # delta_s_topo_2 = np.array(delta_s_topo_2, dtype=np.float64)
        # assert len(delta_s_topo_2) == len(delta_s_topo)
        # delta_s_topo_avg = (delta_s_topo + delta_s_topo_2) / 2
        # delta_s_topo_diff = (delta_s_topo - delta_s_topo_2) / 2
        # if len(delta_s_topo) > 1:
        #     estimate_topo_charge = round((delta_s_topo_avg[0] - (delta_s_topo_avg[1] - 1.5 * delta_s_topo_avg[0]) / 3) / 0.44)
        #     estimate_topo_charge += round((delta_s_topo_diff[0] - (delta_s_topo_diff[1] - 1.5 * delta_s_topo_diff[0]) / 3) / 0.40)
        # elif len(delta_s_topo) > 0:
        #     estimate_topo_charge = round(delta_s_topo_avg[0] / 0.44)
        #     estimate_topo_charge += round(delta_s_topo_diff[0] / 0.40)
        # p_inst["estimate_topo_charge_orig"] = p_inst["estimate_topo_charge"]
        # p_inst["estimate_topo_charge"] = estimate_topo_charge
    # p_inst_list.sort(key=lambda v: v["flow_time"])
    return p_inst_list

@q.timer
def displayln_info_topo_info(topo_info):
    """
    `topo_info` is a `dict`
    """
    info = topo_info
    flow_time = info["flow_time"]
    step_counter = info["step_counter"]
    flow_type = info["flow_type"]
    current_spacing = info["current_spacing"]
    plaq = info["plaq"]
    plaq_min = info["plaq_min"]
    num_instanton = info["num_instanton"]
    topo = info["topo"]
    q.displayln_info(0, f"{flow_type} flow_time={flow_time:9.2f} ; topo={topo:10.3f} ; plaq_min={plaq_min:8.5f} ; plaq={plaq:9.6f} ; num_instanton={num_instanton:11.2f} ; step_counter={step_counter:5} ; current_spacing={current_spacing}")

@q.timer
def displayln_info_p_inst(p_inst):
    """
    e.g.:
    p_inst_list = process_inst_list(inst_list)
    for p_inst in p_inst_list:
        displayln_info_p_inst(p_inst)
    """
    inst_idx = p_inst["inst_idx"]
    inst_raw = p_inst["inst_raw"]
    o_xg_d = p_inst["o_xg_d"]
    o_total_site = p_inst["o_total_site"]
    xg_d = p_inst["xg_d"]
    total_site = p_inst["total_site"]
    current_spacing = p_inst["current_spacing"]
    plaq = p_inst["plaq"]
    flow_time = p_inst["flow_time"]
    flow_time_start = p_inst["flow_time_start"]
    flow_time_end = p_inst["flow_time_end"]
    num_plaq = p_inst["num_plaq"]
    dis_sqr_avg = p_inst["dis_sqr_avg"]
    dis_sqr_max = p_inst["dis_sqr_max"]
    dis_sqr_max_max = p_inst["dis_sqr_max_max"]
    c_dis_sqr_max = p_inst["c_dis_sqr_max"]
    delta_s_topo = p_inst["delta_s_topo"]
    topo_sphere_sum_radius_list = p_inst["topo_sphere_sum_radius_list"]
    estimate_topo_charge = p_inst["estimate_topo_charge"]
    plaq_list = p_inst["plaq_list"]
    s_topo_list = p_inst["s_topo_list"]
    num_plaq_list = p_inst["num_plaq_list"]
    dis_sqr_avg_list = p_inst["dis_sqr_avg_list"]
    dis_sqr_max_list = p_inst["dis_sqr_max_list"]
    c_dis_sqr_list = p_inst["c_dis_sqr_list"]
    closest_sim_inst_info = p_inst["closest_sim_inst_info"]
    closest_inst_info = p_inst["closest_inst_info"]
    idx = inst_idx
    q.displayln_info(0, f"{idx:3}:")
    q.displayln_info(0, f"{idx:3}: o_xg_d={o_xg_d}")
    q.displayln_info(0, f"{idx:3}: o_total_site={o_total_site}")
    q.displayln_info(0, f"{idx:3}: xg_d={xg_d}")
    q.displayln_info(0, f"{idx:3}: total_site={total_site}")
    q.displayln_info(0, f"{idx:3}: current_spacing={current_spacing}")
    q.displayln_info(0, f"{idx:3}: plaq={plaq:.4f}")
    q.displayln_info(0, f"{idx:3}: flow_time={flow_time:.2f}")
    q.displayln_info(0, f"{idx:3}: num_plaq={num_plaq}")
    # q.displayln_info(0, f"{idx:3}: dis_sqr_avg={dis_sqr_avg}")
    q.displayln_info(0, f"{idx:3}: dis_sqr_max={dis_sqr_max}")
    q.displayln_info(0, f"{idx:3}: dis_sqr_max_max={dis_sqr_max_max}")
    q.displayln_info(0, f"{idx:3}: c_dis_sqr_max={c_dis_sqr_max}")
    q.displayln_info(0, f"{idx:3}: flow_time_start={flow_time_start:.2f}")
    q.displayln_info(0, f"{idx:3}: flow_time_end={flow_time_end:.2f}")
    q.displayln_info(0, f"{idx:3}: plaq_list={[ f'{v:.4f}' for v in plaq_list ]}")
    q.displayln_info(0, f"{idx:3}: s_topo_list={[ [ f'{v:.3f}' for v in l ] for l in s_topo_list ]}")
    q.displayln_info(0, f"{idx:3}: topo_sphere_sum_radius_list={topo_sphere_sum_radius_list}")
    q.displayln_info(0, f"{idx:3}: num_plaq_list={num_plaq_list}")
    q.displayln_info(0, f"{idx:3}: dis_sqr_avg_list={[ f'{v:.3f}' for v in dis_sqr_avg_list ]}")
    q.displayln_info(0, f"{idx:3}: dis_sqr_max_list={dis_sqr_max_list}")
    q.displayln_info(0, f"{idx:3}: c_dis_sqr_list={c_dis_sqr_list}")
    q.displayln_info(0, f"{idx:3}: closest_sim_inst_info={closest_sim_inst_info}")
    q.displayln_info(0, f"{idx:3}: closest_inst_info={closest_inst_info}")
    if delta_s_topo is None:
        q.displayln_info(0, f"{idx:3}: delta_s_topo={delta_s_topo}")
    else:
        q.displayln_info(0, f"{idx:3}: delta_s_topo={[ f'{v:.3f}' for v in delta_s_topo ]}")
    q.displayln_info(0, f"{idx:3}: estimate_topo_charge={estimate_topo_charge}")
    # estimate_topo_charge_orig = p_inst["estimate_topo_charge_orig"]
    # if estimate_topo_charge_orig is not None:
    #     q.displayln_info(0, f"{idx:3}: estimate_topo_charge_orig={estimate_topo_charge_orig}")

@q.timer
def displayln_info_p_inst_list(p_inst_list):
    """
    e.g.:
    p_inst_list = process_inst_list(inst_list)
    displayln_info_p_inst_list(p_inst_list)
    """
    tot_topo = 0
    tot_inst = 0
    tot_topo_dict = dict()
    tot_inst_dict = dict()
    for p_inst in p_inst_list:
        displayln_info_p_inst(p_inst)
        inst_idx = p_inst["inst_idx"]
        estimate_topo_charge = p_inst["estimate_topo_charge"]
        current_spacing = p_inst["current_spacing"]
        idx = inst_idx
        tot_topo += estimate_topo_charge
        tot_inst += abs(estimate_topo_charge)
        q.displayln_info(0, f"{idx:3}: tot_inst={tot_inst}")
        q.displayln_info(0, f"{idx:3}: tot_topo={tot_topo}")
        if current_spacing in tot_topo_dict:
            tot_topo_dict[current_spacing] += estimate_topo_charge
            tot_inst_dict[current_spacing] += abs(estimate_topo_charge)
        else:
            tot_topo_dict[current_spacing] = estimate_topo_charge
            tot_inst_dict[current_spacing] = abs(estimate_topo_charge)
    q.displayln_info(0, f"Final:")
    q.displayln_info(0, f"Final: tot_inst_dict={tot_inst_dict}")
    q.displayln_info(0, f"Final: tot_topo_dict={tot_topo_dict}")
    q.displayln_info(0, f"Final: tot_inst={tot_inst}")
    q.displayln_info(0, f"Final: tot_topo={tot_topo}")

def get_tot_topo_count(p_inst_list, flow_time):
    """
    Get the total topological charge by counting the number of instanton minus anti-instanton.
    """
    topo_count = sum([ round(p_inst["estimate_topo_charge"]) for p_inst in p_inst_list if p_inst["flow_time"] > flow_time ])
    return topo_count

def get_tot_inst_count(p_inst_list, flow_time):
    """
    Get the total number of instanton plus anti-instanton.
    """
    inst_count = sum([ abs(round(p_inst["estimate_topo_charge"])) for p_inst in p_inst_list if p_inst["flow_time"] > flow_time ])
    return inst_count

### --------------------------------------

class InstantonMap:

    """
    self.inst_list
    self.active_inst_list
    #
    self.inst_list[inst_idx] = [
        dict(
            flow_time,
            current_spacing,
            total_site,
            o_xg_d,
            o_xg,
            xg_d,
            xg,
            plaq,
            s_topo_list,
            topo_sphere_sum_radius_list,
            num_plaq,
            dis_sqr_avg,
            dis_sqr_max,
            ),
        ...
    ]
    """

    def __init__(self, *, total_site=None, topo_sphere_sum_radius_list=None, dis_sqr_limit=None, threshold=None):
        """
        total_site = q.Coordinate(get_param(job_tag, "total_site"))
        topo_sphere_sum_radius_list = 5
        dis_sqr_limit = 5
        threshold = 0.1
        """
        self.inst_list = list()
        self.active_inst_list = list()
        self.flow_time = 0.0
        self.step_counter = 0
        self.current_spacing = 1
        self.current_flow_time_multiplier = 1
        self.initial_total_site = total_site
        self.total_site = total_site
        self.geo = q.Geometry(self.total_site)
        self.plaq_xg_arr = mk_plaq_xg_arr(self.geo)
        self.origin_coordinate = q.Coordinate()
        # dis_sqr_limit (always in unit of current lattice spacing)
        # limit for the distance square to identify as the same instanton.
        self.dis_sqr_limit = dis_sqr_limit
        # threshold of (1 - plaq) for instanton detection
        # if `plaq < 1 - threshold`, then this plaq will be considered for instanton tunnelling.
        self.threshold = threshold
        self.topo_sphere_sum_radius_list = topo_sphere_sum_radius_list
        self.info_list = []
        self.flow_time_list = []
        self.tot_topo_list = []

    def shift(self, shift):
        """
        If we shift the field with `shift` by
        `field = field.shift(shift)`
        then we also shift `self.origin_coordinate`,
        which represent the current location of the initial origin
        (in unit of the initial lattice space).
        `shift` is in unit of the current lattice space.
        """
        self.origin_coordinate = self.origin_coordinate + shift * self.current_spacing

    def convert_xg(self, xg):
        """
        Convert current `xg` to the coordinate of the initial lattice.
        """
        assert self.initial_total_site is not None
        if isinstance(xg, q.Coordinate):
            o_xg = (xg * self.current_spacing - self.origin_coordinate) % self.initial_total_site
        elif isinstance(xg, q.CoordinateD):
            initial_total_site_d = q.CoordinateD(self.initial_total_site)
            origin_coordinate_d = q.CoordinateD(self.origin_coordinate)
            o_xg = (xg * self.current_spacing - origin_coordinate_d) % initial_total_site_d
        else:
            assert False
        return o_xg

    @q.timer
    def half_lattice(self):
        assert len(self.active_inst_list) == 0
        self.current_spacing *= 2
        self.current_flow_time_multiplier *= 2**4 / 2
        total_site = q.Coordinate([ c // 2 for c in self.total_site ])
        assert total_site * 2 == self.total_site
        assert total_site * self.current_spacing == self.initial_total_site
        self.total_site = total_site
        self.geo = q.Geometry(self.total_site)
        self.plaq_xg_arr = mk_plaq_xg_arr(self.geo)

    def acc_time(self, step_size):
        """
        Call after flow.
        """
        self.flow_time += step_size * self.current_flow_time_multiplier
        self.step_counter += 1

    @q.timer
    def acc_topo_info(self, f_plaq, f_topo=None):
        """
        e.g.
        group_extreme_plaq_xg_list = get_group_extreme_plaq_xg_list(extreme_plaq_xg_list, total_site, dis_sqr_limit)
        f_topo = q.gf_topology_field(gf)
        """
        assert self.total_site is not None
        total_site = self.total_site
        total_site_d = q.CoordinateD(total_site)
        extreme_plaq_xg_list = get_extreme_plaq_xg_list(self.plaq_xg_arr, f_plaq, self.threshold)
        group_extreme_plaq_xg_list = get_group_extreme_plaq_xg_list(extreme_plaq_xg_list, total_site, self.dis_sqr_limit)
        num_group = len(group_extreme_plaq_xg_list)
        if num_group == 0:
            self.inst_list += self.active_inst_list
            self.active_inst_list = []
            return
        for g in group_extreme_plaq_xg_list:
            assert isinstance(g, list)
            assert len(g) > 0
            assert isinstance(g[0], tuple)
            assert len(g[0]) == 3
        plaq_arr = np.array([ g[0][0] for g in group_extreme_plaq_xg_list ], dtype=np.float64)
        xg_d_arr = np.array([ g[0][1] for g in group_extreme_plaq_xg_list ], dtype=np.float64)
        dis_sqr_arr_list = [ np.array([ p[2] for p in g ]) for g in group_extreme_plaq_xg_list ]
        xg_arr = np.array(xg_d_arr, dtype=np.int32)
        if f_topo is not None:
            assert self.topo_sphere_sum_radius_list is not None
            assert isinstance(self.topo_sphere_sum_radius_list, list)
            geo = f_topo.geo
            topo = f_topo.glb_sum()[:].item()
            f_topo_c = q.FieldComplexD(geo, 1)
            f_s_topo = q.FieldRealD(geo, 1)
            f_topo_c[:] = f_topo[:]
            num_radius = len(self.topo_sphere_sum_radius_list)
            s_topo_arr = np.zeros((num_group, num_radius,), dtype=np.float64)
            for idx_radius, radius in enumerate(self.topo_sphere_sum_radius_list):
                f_s_topo_c = q.sphere_sum_field(f_topo_c, radius)
                f_s_topo[:] = f_s_topo_c[:].real
                sp_s_topo = f_s_topo.get_elems_xg(xg_arr)
                s_topo_arr[:, idx_radius] = sp_s_topo[:, 0]
            assert s_topo_arr.shape == (num_group, num_radius,)
        else:
            num_radius = 0
            s_topo_arr = np.array([ None for i in range(num_group) ], dtype=object)
            assert s_topo_arr.shape == (num_group,)
        assert plaq_arr.shape == (num_group,)
        assert xg_d_arr.shape == (num_group, 4,)
        assert xg_arr.shape == (num_group, 4,)
        assert len(dis_sqr_arr_list) == num_group
        prev_active_inst_list = self.active_inst_list
        self.active_inst_list = []
        for idx in range(num_group):
            xg_d = q.CoordinateD(xg_d_arr[idx])
            xg = q.Coordinate(xg_arr[idx])
            o_xg_d = self.convert_xg(xg_d)
            o_xg = self.convert_xg(xg)
            plaq = plaq_arr[idx].item()
            s_topo = s_topo_arr[idx]
            dis_sqr_arr = dis_sqr_arr_list[idx]
            num_plaq = len(dis_sqr_arr)
            assert dis_sqr_arr.shape == (num_plaq,)
            dis_sqr_avg = np.average(dis_sqr_arr).item()
            dis_sqr_max = np.max(dis_sqr_arr).item()
            if s_topo is not None:
                s_topo = s_topo.tolist()
            info_dict = dict()
            info_dict["flow_time"] = self.flow_time
            info_dict["current_spacing"] = self.current_spacing
            info_dict["total_site"] = self.total_site.to_list()
            info_dict["xg_d"] = xg_d.to_list()
            info_dict["xg"] = xg.to_list()
            info_dict["o_xg"] = o_xg.to_list()
            info_dict["o_xg_d"] = o_xg_d.to_list()
            info_dict["plaq"] = plaq
            info_dict["s_topo"] = s_topo
            info_dict["topo_sphere_sum_radius_list"] = self.topo_sphere_sum_radius_list
            info_dict["num_plaq"] = num_plaq
            info_dict["dis_sqr_avg"] = dis_sqr_avg
            info_dict["dis_sqr_max"] = dis_sqr_max
            is_active = False
            for k in range(len(prev_active_inst_list)):
                inst = prev_active_inst_list[k]
                assert len(inst) > 0
                inst_xg_d = q.CoordinateD(inst[-1]["xg_d"])
                dis_sqr = point_d_dis_sqr(inst_xg_d, xg_d, total_site_d)
                if dis_sqr < self.dis_sqr_limit:
                    is_active = True
                    inst.append(info_dict)
                    self.active_inst_list.append(inst)
                    del prev_active_inst_list[k]
                    break
            if not is_active:
                inst = [ info_dict, ]
                self.active_inst_list.append(inst)
        self.inst_list += prev_active_inst_list

### --------------------------------------

@q.timer(is_timer_fork=True)
def compute_inst_map(
        gf,
        *,
        dis_sqr_limit=None,
        threshold=None,
        topo_sphere_sum_radius_list=None,
        plaq_min_threshold=None,
        max_spacing=None,
        flow_step_size_list=None,
        flow_num_step_list=None,
        ):
    """
    e.g.:
    #
    job_tag = "16IH2"
    set_param(job_tag, "traj_list")(list(range(1000, 4000, 10)))
    set_param(job_tag, "load_config_params")(None)
    #
    get_gf = run_gf(job_tag, traj)
    if get_gf is None:
        return
    gf = get_gf()
    #
    inst_map_obj = compute_inst_map(gf) # Note that `gf` should not have boundary be twisted.
    #
    obj = dict()
    obj["inst_list"] = inst_map.inst_list
    obj["info_list"] = inst_map.info_list
    obj["flow_time_list"] = inst_map.flow_time_list
    obj["tot_topo_list"] = inst_map.tot_topo_list
    """
    fname = q.get_fname()
    if topo_sphere_sum_radius_list is None:
        topo_sphere_sum_radius_list = [ 2, 3, 4, 5, 6, ]
    q.displayln_info(0, f"{fname}: topo_sphere_sum_radius_list={topo_sphere_sum_radius_list}")
    if dis_sqr_limit is None:
        # limit to group plaq to identify as instanton
        dis_sqr_limit = 5
    q.displayln_info(0, f"{fname}: dis_sqr_limit={dis_sqr_limit}")
    if threshold is None:
        # threshold for (1 - plaq) to identify as instanton
        threshold = 0.1
    q.displayln_info(0, f"{fname}: threshold={threshold}")
    if plaq_min_threshold is None:
        # ADJUST ME
        # https://inspirehep.net/literature/166110
        plaq_min_threshold = 0.97**(1/16)
        # plaq_min_threshold = 0.99
    q.displayln_info(0, f"{fname}: plaq_min_threshold={plaq_min_threshold}")
    #
    geo = gf.geo
    total_site = geo.total_site
    #
    if max_spacing is None:
        max_spacing = 1
        total_site_arr = total_site.to_numpy()
        while np.all(total_site_arr % (2 * max_spacing) == 0):
            max_spacing *= 2
    #
    if flow_step_size_list is None:
        flow_step_size_list = [ 0.05, 0.1, 0.1, ]
    else:
        assert isinstance(flow_step_size_list, list)
        assert len(flow_step_size_list) == 3
    #
    if flow_num_step_list is None:
        flow_num_step_list = [ 80, 320, 1000000, ]
    else:
        assert isinstance(flow_num_step_list, list)
        assert len(flow_num_step_list) == 3
    #
    inst_map = InstantonMap(
            total_site=total_site,
            topo_sphere_sum_radius_list=topo_sphere_sum_radius_list,
            dis_sqr_limit=dis_sqr_limit,
            threshold=threshold,
            )
    #
    flow_type = None
    plaq = None
    plaq_min = None
    num_instanton = None
    topo = None
    def append_topo_info():
        info = dict()
        info["flow_time"] = inst_map.flow_time
        info["step_counter"] = inst_map.step_counter
        info["flow_type"] = flow_type
        info["current_spacing"] = inst_map.current_spacing
        info["plaq"] = plaq
        info["plaq_min"] = plaq_min
        info["num_instanton"] = num_instanton
        info["topo"] = topo
        inst_map.info_list.append(info)
    #
    flow_type = "Wilson"
    step_size = flow_step_size_list[0]
    num_step = flow_num_step_list[0]
    for i in range(num_step + 1):
        geo = gf.geo
        f_plaq = q.gf_plaq_field(gf)
        plaq_min = np.min(f_plaq.glb_min()[:]).item()
        plaq = f_plaq.glb_sum()[:].sum().item() / geo.total_volume / 6
        num_instanton = ((1 - plaq) * 6 * geo.total_volume) / (8 * np.pi**2 / 6)
        if inst_map.step_counter % 10 == 0:
            f_topo = q.gf_topology_field(gf)
            topo = f_topo.glb_sum()[:].item()
            append_topo_info()
            q.displayln_info(0, f"flow_time={inst_map.flow_time:.3f} ; topo", topo)
        q.displayln_info(0, f"{flow_type} ; {geo.total_site.to_list()} ; i={i} ; step_counter={inst_map.step_counter} ; flow_time={inst_map.flow_time:.3f} ; 1-plaq_min={1-plaq_min:.5f} ; 1-plaq={1-plaq:.5f} ; num_instanton={num_instanton:.2f}")
        if i != num_step:
            gf_flow_topo(gf, step_size)
            inst_map.acc_time(step_size)
    #
    if q.is_test():
        q.json_results_append(f"{fname}: {flow_type} plaq", plaq, 1e-10)
        q.json_results_append(f"{fname}: {flow_type} plaq_min", plaq_min, 1e-8)
        q.json_results_append(f"{fname}: {flow_type} topo", topo, 1e-8)
    # q.timer_display()
    #
    flow_type = "Freeze"
    step_size = flow_step_size_list[1]
    num_step = flow_num_step_list[1]
    for i in range(num_step + 1):
        geo = gf.geo
        f_plaq = q.gf_plaq_field(gf)
        plaq_min = np.min(f_plaq.glb_min()[:]).item()
        plaq = f_plaq.glb_sum()[:].sum().item() / geo.total_volume / 6
        num_instanton = ((1 - plaq) * 6 * geo.total_volume) / (8 * np.pi**2 / 6)
        if inst_map.step_counter % 10 == 0:
            f_topo = q.gf_topology_field(gf)
            topo = f_topo.glb_sum()[:].item()
            append_topo_info()
            q.displayln_info(0, f"flow_time={inst_map.flow_time:.3f} ; topo", topo)
        q.displayln_info(0, f"{flow_type} ; {geo.total_site.to_list()} ; i={i} ; step_counter={inst_map.step_counter} ; flow_time={inst_map.flow_time:.3f} ; 1-plaq_min={1-plaq_min:.5f} ; 1-plaq={1-plaq:.5f} ; num_instanton={num_instanton:.2f}")
        if plaq_min > 0.995:
            break
        if i != num_step:
            gf_flow_topo(gf, step_size, "Freeze")
            inst_map.acc_time(step_size)
    #
    topo = q.gf_topology(gf)
    inst_map.flow_time_list.append(inst_map.flow_time)
    inst_map.tot_topo_list.append(topo)
    q.displayln_info(0, f"flow_time_list={inst_map.flow_time_list} ; tot_topo_list={inst_map.tot_topo_list}")
    if q.is_test():
        q.json_results_append(f"{fname}: {flow_type} plaq", plaq, 1e-10)
        q.json_results_append(f"{fname}: {flow_type} plaq_min", plaq_min, 1e-8)
        q.json_results_append(f"{fname}: {flow_type} topo", topo, 1e-8)
    # q.timer_display()
    #
    flow_type = "Shrink"
    step_size = flow_step_size_list[2]
    num_step = flow_num_step_list[2]
    for i in range(num_step + 1):
        geo = gf.geo
        f_plaq = q.gf_plaq_field(gf)
        plaq_min = np.min(f_plaq.glb_min()[:]).item()
        plaq = f_plaq.glb_sum()[:].sum().item() / geo.total_volume / 6
        num_instanton = ((1 - plaq) * 6 * geo.total_volume) / (8 * np.pi**2 / 6)
        if inst_map.step_counter % 4 == 0:
            f_topo = q.gf_topology_field(gf)
            topo = f_topo.glb_sum()[:].item()
            append_topo_info()
            q.displayln_info(0, f"flow_time={inst_map.flow_time:.3f} ; topo", topo)
        else:
            f_topo = None
        q.displayln_info(0, f"{flow_type} ; {geo.total_site.to_list()} ; i={i} ; step_counter={inst_map.step_counter} ; flow_time={inst_map.flow_time:.3f} ; 1-plaq_min={1-plaq_min:.5f} ; 1-plaq={1-plaq:.5f} ; num_instanton={num_instanton:.2f}")
        inst_map.acc_topo_info(f_plaq, f_topo)
        q.displayln_info(0, f"num_inst={len(inst_map.inst_list)} ; num_active_inst={len(inst_map.active_inst_list)}")
        if num_instanton < 0.1:
            break
        if abs(topo) < 1e-5 and num_instanton < 0.9:
            break
        if plaq_min > plaq_min_threshold and inst_map.current_spacing * 2 <= max_spacing:
            gf_h = q.gf_reduce_half(gf)
            geo = gf.geo
            geo_h = gf_h.geo
            # gf.show_info()
            # gf_h.show_info()
            f_plaq = q.gf_plaq_field(gf)
            plaq_min = np.min(f_plaq.glb_min()[:]).item()
            plaq = f_plaq.glb_sum()[:].sum().item() / geo.total_volume / 6
            num_instanton = ((1 - plaq) * 6 * geo.total_volume) / (8 * np.pi**2 / 6)
            f_plaq_h = q.gf_plaq_field(gf_h)
            plaq_min_h = np.min(f_plaq_h.glb_min()[:]).item()
            plaq_h = f_plaq_h.glb_sum()[:].sum().item() / geo_h.total_volume / 6
            num_instanton_h = ((1 - plaq_h) * 6 * geo_h.total_volume) / (8 * np.pi**2 / 6)
            f_topo = q.gf_topology_field(gf)
            topo = f_topo.glb_sum()[:].item()
            f_topo_h = q.gf_topology_field(gf_h)
            topo_h = f_topo_h.glb_sum()[:].item()
            q.displayln_info(0, f"===================== Half lattice size ===========================")
            q.displayln_info(0, f"flow_time={inst_map.flow_time}")
            q.displayln_info(0, f"min_plaq", plaq_min, plaq_min_h)
            q.displayln_info(0, f"num_instanton", num_instanton, num_instanton_h)
            q.displayln_info(0, f"topo", topo, topo_h)
            q.displayln_info(0, f"===================== Half lattice size ===========================")
            if q.is_test():
                q.json_results_append(f"{fname}: {flow_type} current_spacing={inst_map.current_spacing}")
                q.json_results_append(f"{fname}: {flow_type} plaq", plaq, 1e-8)
                q.json_results_append(f"{fname}: {flow_type} plaq_h", plaq_h, 1e-8)
                q.json_results_append(f"{fname}: {flow_type} plaq_min", plaq_min, 1e-6)
                q.json_results_append(f"{fname}: {flow_type} plaq_min_h", plaq_min_h, 1e-6)
                q.json_results_append(f"{fname}: {flow_type} topo", topo, 1e-4)
                q.json_results_append(f"{fname}: {flow_type} topo_h", topo_h, 1e-4)
            gf = gf_h
            geo = geo_h
            total_site = geo_h.total_site
            inst_map.half_lattice()
            inst_map.flow_time_list.append(inst_map.flow_time)
            inst_map.tot_topo_list.append(topo_h)
            continue
        if i != num_step:
            gf_flow_topo(gf, step_size, "Shrink")
            inst_map.acc_time(step_size)
    obj = dict()
    obj["inst_list"] = inst_map.inst_list
    obj["info_list"] = inst_map.info_list
    obj["flow_time_list"] = inst_map.flow_time_list
    obj["tot_topo_list"] = inst_map.tot_topo_list
    q.displayln_info(0, f"tot_topo_list={inst_map.tot_topo_list}")
    return obj

@q.timer
def displayln_info_inst_map_obj(inst_map_obj):
    """
    e.g.:
    inst_map_obj = compute_inst_map(gf)
    displayln_info_inst_map_obj(inst_map_obj)
    """
    fname = q.get_fname()
    obj = inst_map_obj
    #
    q.displayln_info(0, f"{fname}: {obj.keys()}")
    #
    info_list = obj["info_list"]
    inst_list = obj["inst_list"]
    tot_topo_list = obj["tot_topo_list"]
    flow_time_list = obj["flow_time_list"]
    #
    p_inst_list = process_inst_list(inst_list)
    #
    if q.is_test():
        for p_inst in p_inst_list:
            inst_idx = p_inst["inst_idx"]
            o_xg_d = p_inst["o_xg_d"]
            o_total_site = p_inst["o_total_site"]
            current_spacing = p_inst["current_spacing"]
            plaq = p_inst["plaq"]
            flow_time = p_inst["flow_time"]
            delta_s_topo = p_inst["delta_s_topo"]
            estimate_topo_charge = p_inst["estimate_topo_charge"]
            q.json_results_append(f"inst_idx={inst_idx} ; o_xg_d={o_xg_d} ; o_total_site={o_total_site} ; current_spacing={current_spacing}")
            q.json_results_append(f"plaq", plaq, 1e-6)
            q.json_results_append(f"flow_time", flow_time, 1e-3)
            q.json_results_append(f"delta_s_topo", np.array(delta_s_topo, dtype=np.float64), 1e-6)
            q.json_results_append(f"estimate_topo_charge", estimate_topo_charge, 1e-4)
    #
    s_p_inst_list = sorted(p_inst_list, key=lambda v: -v["flow_time"])
    #
    q.displayln_info(0, f"")
    for idx, (flow_time, tot_topo,) in enumerate(zip(flow_time_list, tot_topo_list)):
        topo_count = get_tot_topo_count(p_inst_list, flow_time)
        inst_count = get_tot_inst_count(p_inst_list, flow_time)
        current_spacing = 2**idx
        q.displayln_info(0, f"{idx:3}: flow_time={flow_time:9.2f} ; tot_topo={tot_topo:10.3f} ; topo_count={topo_count:4} ; inst_count={inst_count:4} ; current_spacing={current_spacing}")
    #
    q.displayln_info(0, f"")
    #
    for topo_info in info_list:
        flow_type = topo_info["flow_type"]
        flow_time = topo_info["flow_time"]
        while len(s_p_inst_list) > 0:
            if s_p_inst_list[-1]["flow_time"] <= flow_time:
                displayln_info_p_inst(s_p_inst_list.pop())
            else:
                break
        displayln_info_topo_info(topo_info)
        topo_count = get_tot_topo_count(p_inst_list, flow_time)
        inst_count = get_tot_inst_count(p_inst_list, flow_time)
        q.displayln_info(0, f"{flow_type} flow_time={flow_time:9.2f} ; topo_count={topo_count:4} ; inst_count={inst_count:4}")
    #
    q.displayln_info(0, f"")
    #
    tot_topo_direct = round(tot_topo_list[0])
    tot_topo_count = get_tot_topo_count(p_inst_list, 0)
    tot_inst_count = get_tot_inst_count(p_inst_list, 0)
    #
    q.displayln_info(0, f"{fname}: tot_topo_direct={tot_topo_direct} ; tot_topo_count={tot_topo_count}")
    if tot_topo_direct == tot_topo_count:
        q.displayln_info(0, f"{fname}: Topological charge matched.")
    else:
        q.displayln_info(-1, f"{fname}: WARNING: Topological charge does not match.")

### --------------------------------------
