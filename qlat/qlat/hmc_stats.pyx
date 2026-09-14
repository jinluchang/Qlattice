# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

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

from qlat_utils.all cimport *
from . cimport everything as cc
from .hmc cimport GaugeMomentum
from .qcd cimport GaugeField

import qlat_utils as q
import cqlat as c

def get_gm_force_magnitudes(gm_force, n_elems):
    return c.get_gm_force_magnitudes(gm_force, n_elems)

def display_gm_force_magnitudes(GaugeMomentum gm_force, cc.Int n_elems):
    cc.display_gm_force_magnitudes(gm_force.xxx().val(), n_elems)

def save_gm_force_magnitudes_list(fn):
    q.mk_file_dirs_info(fn)
    cc.save_gm_force_magnitudes_list(fn)

def display_gauge_field_info_table_with_wilson_flow(
    fn_gf_info, fn_wilson_flow_energy, GaugeField gf, cc.RealD flow_time,
    cc.Int flow_steps, cc.Int steps, cc.RealD c1=0.0
):
    cc.display_gauge_field_info_table_with_wilson_flow(
        fn_gf_info, fn_wilson_flow_energy, gf.xxx().val(),
        flow_time, flow_steps, steps, c1
    )
