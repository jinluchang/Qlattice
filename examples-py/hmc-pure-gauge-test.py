#!/usr/bin/env python3

import qlat as q

@q.timer_verbose
def test_hmc(total_site, ga):
    #
    q.qmkdir_info("results")
    q.qmkdir_info("results/gf_info")
    q.qmkdir_info("results/wilson_flow_energy_info")
    #
    geo = q.Geometry(total_site)
    rs = q.RngState(f"test_hmc-{total_site.to_list()}")
    gf = q.GaugeField(geo)
    q.set_unit(gf)
    traj = 0
    for i in range(4):
        traj += 1
        q.run_hmc_pure_gauge(
            gf, ga, traj, rs.split("run_hmc_pure_gauge"), is_always_accept=True
        )
        plaq_avg = q.gf_avg_plaq(gf)
        geo.total_volume * 6.0 * (1.0 - plaq_avg)
        wilson_loop = q.gf_avg_wilson_loop_normalized_tr(gf, 1, 1)
        q.json_results_append(
            f"beta={ga.beta()}_c1={ga.c1()} traj={traj} plaq", plaq_avg, 1e-10
        )
        q.json_results_append(
            f"beta={ga.beta()}_c1={ga.c1()} traj={traj} wilson_loop", wilson_loop, 1e-10
        )
        if traj % 2 == 0:
            q.display_gauge_field_info_table_with_wilson_flow(
                f"results/gf_info/traj={traj}.lat",
                f"results/wilson_flow_energy_info/traj={traj}.lat",
                gf,
                0.1,
                5,
                2,
            )

@q.timer_verbose
def main():
    total_site = q.Coordinate(
        [
            4,
            4,
            4,
            8,
        ]
    )
    ga = q.GaugeAction(2.13, -0.331)
    test_hmc(total_site, ga)
    ga = q.GaugeAction(5.5, 0.0)
    test_hmc(total_site, ga)

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 2, 2, 2],
    [2, 2, 2, 2],
    [2, 2, 2, 4],
]

q.begin_with_mpi(size_node_list)

# q.show_machine()

q.qremove_all_info("results")

main()

q.timer_display()

q.check_log_json(__file__)
q.end_with_mpi()

q.displayln_info("CHECK: finished successfully.")
