#!/usr/bin/env python3

import numpy as np
import qlat as q
import qlat_cps as qc
import qlat_gpt as qg
import gpt as g
import qlat_gpt_cps as qgc

from qlat_scripts.v1 import *
from sparse_prop_selection import *
from gen_data_cps_fake import *
from sparse_wsrc_prop_save_from_cps import *
from set_params import *

job_tag = "test-8nt16"
# job_tag = "96I"

total_site = q.Coordinate(get_param(job_tag, "total_site"))

qgc.begin_with_gpt_cps(total_site)

gen_all_data_cps_fake(job_tag)

for traj in get_param(job_tag, "traj_list"):
    run_job(job_tag, traj)

q.timer_display()

qgc.end_with_gpt_cps()

q.displayln_info(f"CHECK: finished successfully.")
