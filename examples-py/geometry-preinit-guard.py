#!/usr/bin/env python3

# Pre-initialisation guard for Geometry(total_site).
#
# ``q.Geometry(q.Coordinate(total_site))`` divides total_site by the global
# size_node, which is unknown until q.begin_with_mpi() (or q.begin()) has run.
# It used to die with SIGFPE (exit 136); it must now raise RuntimeError.
#
# The child interpreter is used because the failure mode was a native signal,
# which would kill this test process instead of being reported.

import os
import subprocess
import sys

import numpy as np

import qlat as q

child_code = (
    "import sys\n"
    "import qlat as q\n"
    "try:\n"
    "    q.Geometry(q.Coordinate([4, 4, 4, 4]))\n"
    "except RuntimeError:\n"
    "    print('CHILD: RuntimeError')\n"
    "except BaseException as e:\n"
    "    print('CHILD: OTHER', type(e).__name__, e)\n"
    "else:\n"
    "    print('CHILD: NO-EXCEPTION')\n"
    "print('CHILD-DONE', sys.version.split()[0])\n"
)

child_env = dict(os.environ)
# inherit the parent's import path (PYTHONPATH from the environment alone is not
# enough once we are running from a copied work directory)
child_env["PYTHONPATH"] = os.pathsep.join(
    p for p in sys.path if p and os.path.isdir(p)
)
r = subprocess.run(
    [sys.executable, "-c", child_code],
    cwd=os.path.dirname(os.path.abspath(__file__)),
    env=child_env,
    capture_output=True,
    text=True,
    timeout=300,
)
child_outcome = 9.0
child_ok = False
for line in r.stdout.splitlines():
    if line.startswith("CHILD: "):
        word = line[7:]
        if word == "RuntimeError":
            child_outcome = 1.0
            child_ok = True
        elif word == "NO-EXCEPTION":
            child_outcome = 0.0
        else:
            child_outcome = 9.0
q.json_results_append("pre-init Geometry(total_site) outcome", child_outcome)
q.json_results_append("pre-init Geometry(total_site) ok", float(child_ok))
q.json_results_append("pre-init child exit code", float(r.returncode))
if child_outcome != 1.0:
    print("CHILD-STDOUT:", r.stdout)
    print("CHILD-STDERR:", r.stderr)

# the guard must not fire once qlat is initialised
q.begin_with_mpi([q.Coordinate([1, 1, 1, 1])])
geo = q.Geometry(q.Coordinate([4, 4, 4, 4]))
q.json_results_append("initialised Geometry(total_site)", np.array(geo.total_site.to_list(), dtype=float))

# Geometry() stays usable without initialisation, and so does the explicit
# (id_node, size_node, node_site) form, which does not consult the global geon
geo0 = q.Geometry()
q.json_results_append("Geometry() total_volume", float(geo0.total_volume))
geo1 = q.Geometry(0, q.Coordinate([1, 1, 1, 1]), q.Coordinate([4, 4, 4, 4]))
q.json_results_append("explicit Geometry(id, size, node)", np.array(geo1.total_site.to_list(), dtype=float))

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-14)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
