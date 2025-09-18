r"""

# Flow and lattice scale determination

Generalize standard Wilson flow to spatial flow.

Also, we can define additional flow observable ($W_2(t)$).

$$
\ba
W_0(t) =& t^2 \la E(t) \ra
\\
W_1(t) =& W(t) = t \frac{d}{dt} (t^2 \la E(t) \ra)
\\
W_2(t) =& -t^3 \frac{d}{dt}\la E(t) \ra
\ea
$$

https://arxiv.org/abs/1006.4518 for $W_0(t)$.
https://arxiv.org/abs/1203.4469 for $W_1(t)$.
"""

### --------------------------------------

import numpy as np

class q:

    from qlat_utils import (
            timer,
            )

class q:

    from qlat_utils import (
            timer,
            displayln_info,
            TimerFork,
            timer_display,
            set_verbose_level,
            cache_call,
            get_fname,
            does_file_exist_qar_sync_node,
            json_results_append,
            save_pickle_obj,
            RngState,
            is_test,
            get_data_sig_arr,
            )

    from .c import (
            FieldRealD,
            GaugeField,
            gf_plaq_flow_force,
            gf_plaq_field,
            gf_energy_density_dir_field,
            gf_evolve,
            )

### --------------------------------------

is_spatial_default = False
t_dir_default = 3
integrator_type_default = "runge-kutta"
flow_time_default = 0.0

@q.cache_call(maxsize=8)
@q.timer
def get_plaq_factor_for_gf_scale_flow(geo, is_spatial, t_dir):
    """
    return plaq_factor
    `plaq_factor` can be used in `gf_plaq_flow_force`.
    ```
    plaq_factor[:].shape == (local_volume, 6,)
    ```
    The 6 elements of plaq_factor represent plaq in directions:
    0:xy, 1:xz, 2:xt, 3:yz, 4:yt, 5:zt
    """
    plaq_factor = q.FieldRealD(geo, 6)
    if is_spatial:
        plaq_factor.set_zero()
        pfv = plaq_factor[:]
        if t_dir == 3:
            # t dir
            pfv[:, 0] = 1
            pfv[:, 1] = 1
            pfv[:, 3] = 1
        elif t_dir == 0:
            # z dir
            pfv[:, 3] = 1
            pfv[:, 4] = 1
            pfv[:, 5] = 1
        elif t_dir == 1:
            # y dir
            pfv[:, 1] = 1
            pfv[:, 2] = 1
            pfv[:, 5] = 1
        elif t_dir == 2:
            # z dir
            pfv[:, 0] = 1
            pfv[:, 2] = 1
            pfv[:, 4] = 1
        else:
            assert False
    else:
        plaq_factor.set_unit()
    return plaq_factor

@q.timer
def gf_flow_scale(gf, step_size, *, is_spatial=None, t_dir=None, integrator_type=None):
    """
    Modify `gf` in place.
    Default Wilson flow with Runge-Kutta integrator.
    ```
    is_spatial in [ False, True, ]
    t_dir in [ 0, 1, 2, 3, ]
    integrator_type in [ "euler", "runge-kutta", ]
    ```
    Runge-Kutta scheme: http://arxiv.org/abs/1006.4518v3
    """
    fname = q.get_fname()
    if is_spatial is None:
        is_spatial = is_spatial_default
    if t_dir is None:
        t_dir = t_dir_default
    if integrator_type is None:
        integrator_type = integrator_type_default
    geo = gf.geo
    plaq_factor = get_plaq_factor_for_gf_scale_flow(geo, is_spatial, t_dir)
    if integrator_type == "euler":
        gm_force = q.gf_plaq_flow_force(gf, plaq_factor)
        q.gf_evolve(gf, gm_force, step_size)
    elif integrator_type == "runge-kutta":
        z = q.gf_plaq_flow_force(gf, plaq_factor)
        z *= 1.0 / 4.0
        q.gf_evolve(gf, z, step_size)
        z *= 17.0 / 9.0
        zp = q.gf_plaq_flow_force(gf, plaq_factor)
        zp *= 8.0 / 9.0
        zp -= z
        q.gf_evolve(gf, zp, step_size)
        z = q.gf_plaq_flow_force(gf, plaq_factor)
        z *= 3.0 / 4.0
        z -= zp
        q.gf_evolve(gf, z, step_size)
    else:
        raise Exception(f"{fname}: integrator_type={integrator_type}")

@q.timer
def gf_plaq_tslice(gf, *, t_dir=None):
    """
    return plaq_arr
    ```
    t_size == geo.total_site[t_dir]
    plaq_arr.shape == (t_size, 6,)
    ```
    """
    if t_dir is None:
        t_dir = t_dir_default
    geo = gf.geo
    total_volume = geo.total_volume
    t_size = geo.total_site[t_dir]
    f_plaq = q.gf_plaq_field(gf)
    plaq_arr = f_plaq.glb_sum_tslice(t_dir=t_dir)[:] / (total_volume / t_size)
    plaq_arr = np.array(plaq_arr, dtype=np.float64)
    return plaq_arr

@q.timer
def gf_energy_density_dir_tslice(gf, *, t_dir=None):
    r"""
    return energy_density_dir_arr
    ```
    t_size == geo.total_site[t_dir]
    energy_density_dir_arr.shape == (t_size, 6,)
    ```
    Similar to `gf_plaq_tslice`
    `gf_energy_density_dir_tslice(gf)[:]` \approx `6 * (1 - gf_plaq_tslice(gf)[:])`
    """
    if t_dir is None:
        t_dir = t_dir_default
    geo = gf.geo
    total_volume = geo.total_volume
    t_size = geo.total_site[t_dir]
    f_edd = q.gf_energy_density_dir_field(gf)
    energy_density_dir_arr = f_edd.glb_sum_tslice(t_dir=t_dir)[:] / (total_volume / t_size)
    energy_density_dir_arr = np.array(energy_density_dir_arr, dtype=np.float64)
    return energy_density_dir_arr

@q.timer(is_timer_fork=True)
def gf_flow_record(
        gf, step_size, num_step,
        *,
        flow_time=None,
        is_spatial=None,
        t_dir=None,
        integrator_type=None,
        ):
    """
    return obj_record 
    obj_record = dict(
        info_list=[
            dict(flow_time, plaq_tslice, energy_density_dir_tslice,),
            ...
        ],
        params=dict(step_size, num_step, flow_time, is_spatial, t_dir, integrator_type,),
    )
    Modify `gf` in place.
    Default Wilson flow with Runge-Kutta integrator.
    ```
    is_spatial in [ False, True, ]
    t_dir in [ 0, 1, 2, 3, ]
    integrator_type in [ "euler", "runge-kutta", ]
    ```
    """
    if flow_time is None:
        flow_time = flow_time_default
    params = dict(
            step_size=step_size,
            num_step=num_step,
            flow_time=flow_time,
            is_spatial=is_spatial,
            t_dir=t_dir,
            integrator_type=integrator_type,
            )
    info_list = []
    for i in range(num_step):
        q.gf_flow_scale(
                gf, step_size,
                is_spatial=is_spatial, t_dir=t_dir,
                integrator_type=integrator_type,
                )
        gf.unitarize()
        flow_time += step_size
        plaq_tslice = gf_plaq_tslice(gf)
        energy_density_dir_tslice = gf_energy_density_dir_tslice(gf)
        info = dict(
            flow_time=flow_time,
            plaq_tslice=plaq_tslice,
            energy_density_dir_tslice=energy_density_dir_tslice,
        )
        info_list.append(info)
    obj_record = dict(
            info_list=info_list,
            params=params,
            )
    return obj_record

default_run_flow_scale_params = dict(
        step_size=0.05,
        num_step=400,
        is_spatial=False,
        t_dir=3,
        integrator_type="runge-kutta",
        )

@q.timer(is_timer_fork=True)
def run_flow_scale(fn_out, *, get_gf=None, fn_gf=None, params=None):
    fname = q.get_fname()
    if not fn_out.endswith(".pickle"):
        q.displayln_info(-1, f"{fname}: WARNING: '{fn_out}' does not endswith '.pickle'. Skip this file.")
        return
    if q.does_file_exist_qar_sync_node(fn_out):
        q.displayln_info(-1, f"{fname}: WARNING: '{fn_out}' for '{fn_gf}' already exist.")
        return
    if get_gf is None:
        if not q.does_file_exist_qar_sync_node(fn_gf):
            q.displayln_info(-1, f"{fname}: WARNING: '{fn_gf}' does not exist. Skip this file.")
            return
    else:
        assert fn_gf is None
    if params is None:
        params = q.default_run_flow_scale_params.copy()
    else:
        params = q.default_run_flow_scale_params | params
    q.json_results_append(f"{fname}: Start compute flow scale info fn='{fn_out}' for '{fn_gf}'")
    if get_gf is None:
        gf = q.GaugeField()
        gf.load(fn_gf)
    else:
        gf = get_gf().copy()
    step_size = params["step_size"]
    num_step = params["num_step"]
    is_spatial = params["is_spatial"]
    t_dir = params["t_dir"]
    integrator_type = params["integrator_type"]
    obj_record = q.gf_flow_record(
            gf, step_size, num_step,
            is_spatial=is_spatial, t_dir=t_dir, integrator_type=integrator_type)
    q.save_pickle_obj(obj_record, fn_out)
    if q.is_test():
        q.json_results_append(f"{fname}: params={obj_record['params']}")
        q.json_results_append(f"{fname}: flow_time", obj_record['info_list'][-1]['flow_time'])
        rs = q.RngState()
        q.json_results_append(
                f"{fname}: sig(plaq_tslice)",
                q.get_data_sig_arr(obj_record['info_list'][-1]['plaq_tslice'], rs, 3))
        q.json_results_append(
                f"{fname}: sig(energy_density_dir_tslice)",
                q.get_data_sig_arr(obj_record['info_list'][-1]['energy_density_dir_tslice'], rs, 3))
    q.json_results_append(f"{fname}: End compute flow scale info fn='{fn_out}' for '{fn_gf}'")

### --------------------------------------

(
        q.get_plaq_factor_for_gf_scale_flow,
        q.gf_flow_scale,
        q.gf_plaq_tslice,
        q.gf_energy_density_dir_tslice,
        q.gf_flow_record,
        q.run_flow_scale,
        q.default_run_flow_scale_params,
        ) = (
                get_plaq_factor_for_gf_scale_flow,
                gf_flow_scale,
                gf_plaq_tslice,
                gf_energy_density_dir_tslice,
                gf_flow_record,
                run_flow_scale,
                default_run_flow_scale_params,
                )
