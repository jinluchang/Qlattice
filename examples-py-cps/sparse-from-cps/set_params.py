import numpy as np
import qlat as q
import qlat_gpt as qg
import gpt as g

from qlat_scripts.v1 import *
from sparse_prop_selection import *

load_path_list[:] = [
    "results",
    "results-fake",
]

set_param_field_selection_rate("96I", 1 / 64, 4096)
set_param_field_selection_rate("test-8nt16", 1 / 16, 128)

set_param("96I", "traj_list")([ traj for traj in range(500, 2000, 10) ])

set_param("test-8nt16", "traj_list")([ traj for traj in range(1000, 2000, 100) ])