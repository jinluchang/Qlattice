"""
Module ``qlat.hmc_stats``
==========================\n
HMC diagnostic and analysis utilities.\n
Provides helpers for inspecting molecular-dynamics force magnitudes during
HMC trajectories and for generating gauge-field information tables with
Wilson-flow analysis.\n
Documentation: ``docs/qlat/qlat_hmc_stats.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat.hmc import *

import qlat.c as c
import qlat_utils as q

def get_gm_force_magnitudes(gm_force, n_elems):
    return c.get_gm_force_magnitudes(gm_force, n_elems)

def display_gm_force_magnitudes(gm_force, n_elems):
    c.display_gm_force_magnitudes(gm_force, n_elems)

def save_gm_force_magnitudes_list(fn):
    q.mk_file_dirs_info(fn)
    c.save_gm_force_magnitudes_list(fn)

def display_gauge_field_info_table_with_wilson_flow(
    fn_gf_info, fn_wilson_flow_energy, gf, flow_time, flow_steps, steps, c1=0.0
):
    c.display_gauge_field_info_table_with_wilson_flow(
        fn_gf_info, fn_wilson_flow_energy, gf, flow_time, flow_steps, steps, c1
    )
