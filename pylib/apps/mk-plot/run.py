#!/usr/bin/env python3

import qlat.qplot as q
# import qlat as q
import numpy as np
import os

q.qremove_all("plots")

x = np.arange(31) * (6 / 30) - 3

y = np.sin(x)

yerr = 0.1 / (1 + x**2)

y1 = np.cos(x)

y1err = 0.1 / (1 + x**2)

dt = q.azip(x, y, yerr)

dt1 = q.azip(x, y1, y1err)

dt2 = q.read_datatable(q.show_datatable(dt))

assert np.linalg.norm(dt2 - dt) < 1e-8

q.save_datatable(dt, "plots/table.txt")

fn_img = q.plot_save(
        fn = "plots/example.png",
        dts = {
            "table.txt" : dt,
            "table-1.txt" : dt1,
            },
        cmds = [
            "set title 'Example'",
            "set key tm",
            "set size 0.6, 1.0",
            "set xlabel '$x$'",
            "set ylabel '$y$'",
            ],
        lines = [
            "plot [:] [:]",
            "0 lt 0 not",
            "'table.txt' lt 1 w yerrorb t '$\\sin$'",
            "'table-1.txt' lt 2 w yerrorb t '$\\cos$'",
            ],
        # is_run_make = False, # default is True
        # is_display = True , # default is False
        )

fn = os.path.join("plots/example.pyplot.dir", "table.txt")

dt3 = q.load_datatable(fn)

assert np.linalg.norm(dt3 - dt) < 1e-8
