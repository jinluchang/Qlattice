import argparse
from pprint import pformat

import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 3],
    [1, 1, 1, 4],
    [1, 1, 1, 6],
    [1, 1, 1, 8],
    [1, 2, 2, 4],
    [2, 2, 2, 4],
    [2, 2, 2, 4],
]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Topological charge measurement based on Wilson flow with Qlattice"
    )
    parser.add_argument(
        "--source",
        type=str,
        default=None,
        help="Source gauge field config file",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output path for results",
    )
    parser.add_argument(
        "--density_field",
        action="store_true",
        default=False,
        help="Save density field",
    )
    parser.add_argument(
        "--flow_time",
        type=int,
        default=6,
        help="Flow time (default: 6)",
    )
    parser.add_argument(
        "--flow_n_step_per_unit_time",
        type=int,
        default=80,
        help="Flow steps per unit time (default: 80)",
    )
    args, _ = parser.parse_known_args()
    return args

if __name__ == "__main__":
    sys_args = parse_args()

    q.begin_with_mpi(size_node_list)

    q.displayln_info(
        "Topological charge measurement based on Wilson flow with Qlattice"
    )
    q.displayln_info("by Luchang Jin")
    q.displayln_info("2024/05/20")

    p_source = sys_args.source
    p_output = sys_args.output
    flow_time = sys_args.flow_time
    flow_n_step = sys_args.flow_n_step_per_unit_time
    is_density_field = sys_args.density_field

    info_path = p_output

    if is_density_field:
        assert p_output is not None
        density_field_path = info_path
    else:
        density_field_path = None

    smear_info_list = [
        [
            1.0 / flow_n_step,
            flow_n_step,
            0.0,
            "runge-kutta",
        ],
    ] * flow_time

    energy_derivative_info = [
        1.0 / flow_n_step,
        0.0,
        "runge-kutta",
    ]

    def load():
        if p_source is None:
            q.displayln_info(
                "Need to provide source file with '--source filename'. Use a sample gauge field for now."
            )
            total_site = q.Coordinate(
                [
                    4,
                    4,
                    4,
                    8,
                ]
            )
            geo = q.Geometry(total_site)
            gf = q.GaugeField(geo)
            rs = q.RngState("seed")
            gf.set_rand(rs.split("gf-init"), 0.5, 10)
        else:
            gf = q.GaugeField()
            gf.load(p_source)
        q.clear_all_caches()
        q.clear_mem_cache()
        return gf

    gf = load()

    (
        topo_list,
        energy_list,
    ) = q.smear_measure_topo(
        gf,
        smear_info_list=smear_info_list,
        energy_derivative_info=energy_derivative_info,
        info_path=info_path,
        density_field_path=density_field_path,
    )

    if info_path is None:
        q.displayln_info(
            "To save the result, use '--output path'. Print to screen for now."
        )
        q.displayln_info(pformat(topo_list))
        q.displayln_info(pformat(energy_list))

    q.timer_display()

    q.end_with_mpi()
