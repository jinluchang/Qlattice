#!/usr/bin/env python3

import qlat.qplot as q
import numpy as np
import os

q.remove_folder("plots")

x = np.arange(30) * (6 / 30) - 3

y = np.sin(x)

err = 1 / (1 + x**2)

dt = q.azip(x, y, err)

dt2 = q.read_datatable(q.show_datatable(dt))

assert np.linalg.norm(dt2 - dt) < 1e-8

q.save_datatable(dt, "plots/table.txt")

fn_img = q.plot_save(
        fn = "plots/example.png",
        dts = { "table.txt" : dt },
        lines = [
            "plot []",
            "x not",
            "'table.txt' w yerrorb not",
            ],
        # is_display = True,
        )

fn = os.path.join("plots/example.pyplot.dir", "table.txt")

dt3 = q.load_datatable(fn)

assert np.linalg.norm(dt3 - dt) < 1e-8
