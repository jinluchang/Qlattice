#!/usr/bin/env python3

import numpy as np
import qlat as q

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")

q.qplot.plot_save_display_width = 500

x = np.arange(31) * (6 / 30) - 3
y = np.cos(x)
yerr = 0.1 / (1 + x**2)
dts = {
    "table.txt": q.azip(x, y, yerr),
}

if q.get_id_node() == 0:
    q.plot_save(
        fn = "results/plot.png",
        dts = dts,
        cmds = [
            "set size 0.8, 1.0",
            "set key tm",
            "set xlabel '$x$'",
            "set ylabel '$y$'",
        ],
        lines = [
            "plot [-3:3] [-1.5:1.5]",
            "0 not",
            "sin(x) w l t '$y = \\sin(x)$'",
            "'table.txt' w yerrorb t '$y = \\cos(x)$'",
        ],
    )

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
