import numpy as np
import qlat as q
import qlat_cps as qc
import qlat_gpt as qg
import gpt as g

from qlat_scripts.v1 import *

@q.timer_verbose
def gen_gauge_transform_cps_fake(job_tag, traj):
    fname = q.get_fname()
    path = f"{job_tag}/gauge-transform-cps/traj-{traj}.gfix"
    if get_load_path(path) is not None:
        q.displayln_info(-1, f"{fname}: {job_tag}/{traj} already done.")
        return
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    rs = q.RngState(f"{fname}/{job_tag}/{traj}")
    gt = q.GaugeTransform(geo)
    gt.set_rand(rs)
    gt.unitarize()
    gt.save_cps(f"results-fake/{path}")
    q.displayln_info(-1, f"{fname}: {job_tag}/{traj} finish.")

@q.timer_verbose
def gen_prop_wsrc_cps_fake(job_tag, traj, inv_type):
    """
    Generate some fake data.
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", "charm", ]
    inv_type_name = inv_type_name_list[inv_type]
    path = f"{job_tag}/prop-wsrc-full-cps-{inv_type_name}/traj-{traj}"
    if get_load_path(f"{path}/checkpoint.txt"):
        q.displayln_info(-1, f"{fname}: {job_tag}/{traj}/{inv_type_name} already done.")
        return
    with q.TimerFork(max_call_times_for_always_show_info=0, verbose=1, show_display=True):
        total_site = q.Coordinate(get_param(job_tag, "total_site"))
        geo = q.Geometry(total_site)
        num_prop_sel = total_site[3] // 2
        num_prop_exact_sel = 2
        rs = q.RngState(f"{fname}-sel/{job_tag}/{traj}")
        tslice_list = []
        for i in range(num_prop_sel):
            tslice = rs.rand_gen() % total_site[3]
            if tslice in tslice_list:
                continue
            tslice_list.append(tslice)
        tslice_exact_list = []
        for i in range(num_prop_exact_sel):
            tslice = rs.select(tslice_list)
            if tslice in tslice_exact_list:
                continue
            tslice_exact_list.append(tslice)
        rs = q.RngState(f"{fname}-prop/{job_tag}/{traj}")
        prop = q.Prop(geo)
        for tslice in tslice_exact_list:
            inv_acc = 2
            tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
            prop.set_rand_g(rs.split(f"{tslice}-exact"), 0.0, 1.0)
            qc.save_cps_prop_double(
                prop,
                f"results-fake/{path}/{tag}.cps-prop",
            )
        for tslice in tslice_list:
            inv_acc = 1
            tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
            prop.set_rand_g(rs.split(f"{tslice}-sloppy"), 0.0, 1.0)
            qc.save_cps_prop_double(
                prop,
                f"results-fake/{path}/{tag}.cps-prop",
            )
        q.qtouch_info(f"results-fake/{path}/checkpoint.txt")
        q.displayln_info(-1, f"{fname}: {job_tag}/{traj}/{inv_type_name} finish.")

@q.timer_verbose
def gen_all_data_cps_fake(job_tag):
    with q.TimerFork(max_call_times_for_always_show_info=0, verbose=1, show_display=True):
        fname = q.get_fname()
        num_traj = 2
        rs = q.RngState(f"{fname}-traj-sel/{job_tag}")
        for i in range(num_traj):
            traj = rs.select(get_param(job_tag, "trajs"))
            gen_gauge_transform_cps_fake(job_tag, traj)
            for inv_type in [ 0, 1, ]:
                gen_prop_wsrc_cps_fake(job_tag, traj, inv_type)
