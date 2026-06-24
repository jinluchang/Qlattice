#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g

from qlat_scripts.v1.rbc_ukqcd import (
    get_param_fermion,
    get_ls_from_fermion_params,
    get_param_lanc,
    get_param_clanc,
    get_param_cg_mp_maxiter,
    mk_pc_parity,
    mk_pc_ne,
    mk_quark_matrix,
    mk_gpt_inverter,
    mk_inverter,
    get_inv,
)
from qlat_scripts.v1 import get_param

qg.begin_with_gpt()

job_tag = "test-4nt8"
inv_type = 0
inv_acc = 0

# get_ls_from_fermion_params
fp_Ls = {"Ls": 12, "M5": 1.8}
q.json_results_append(
    "get_ls_from_fermion_params(Ls)", get_ls_from_fermion_params(fp_Ls)
)
fp_omega = {"omega": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]}
q.json_results_append(
    "get_ls_from_fermion_params(omega)", get_ls_from_fermion_params(fp_omega)
)

# get_param_fermion
fp = get_param_fermion(job_tag, inv_type, inv_acc)
assert fp is not None
assert fp["M5"] == 1.8
q.json_results_append("get_param_fermion mass", fp["mass"])
q.json_results_append("get_param_fermion M5", fp["M5"])
fp_fallback = get_param_fermion(job_tag, inv_type, 99)
assert fp_fallback is not None
q.json_results_append("get_param_fermion fallback mass", fp_fallback["mass"])

# get_param_lanc
lanc = get_param_lanc(job_tag, inv_type, inv_acc)
assert lanc is not None
assert "fermion_params" in lanc
assert "irl_params" in lanc
assert lanc["fermion_params"] == fp
q.json_results_append("get_param_lanc Nstop", lanc["irl_params"]["Nstop"])

# get_param_clanc
clanc = get_param_clanc(job_tag, inv_type, inv_acc)
assert clanc is not None
assert clanc["nbasis"] <= lanc["irl_params"]["Nstop"]
q.json_results_append("get_param_clanc nbasis", clanc["nbasis"])

# get_param_cg_mp_maxiter
maxiter = get_param_cg_mp_maxiter(job_tag, inv_type, inv_acc)
assert maxiter > 0
q.json_results_append("get_param_cg_mp_maxiter", maxiter)

# GPT fermion setup
total_site = q.Coordinate(get_param(job_tag, "total_site"))
geo = q.Geometry(total_site)
rs = q.RngState("seed")
gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.05, 2)
q.json_results_append("gf.plaq()", gf.plaq(), 1e-12)

parity = mk_pc_parity(job_tag, inv_type, inv_acc)
assert parity == g.odd

gpt_gf = qg.gpt_from_qlat(gf)
qm = mk_quark_matrix(job_tag, gpt_gf, inv_type, inv_acc)
assert qm is not None

pc_ne = mk_pc_ne(job_tag, inv_type, inv_acc)
assert pc_ne is not None

# mk_gpt_inverter
inv = mk_gpt_inverter(gf, job_tag, inv_type, inv_acc, n_grouped=q.get_num_node())
assert inv is not None

src = q.Prop(geo)
src.set_rand(rs.split("prop-src"))
q.json_results_append("src.qnorm()", src.qnorm(), 1e-10)
sol = inv * src
q.json_results_append("sol.qnorm()", sol.qnorm(), 1e-10)

# mk_inverter alias
inv2 = mk_inverter(gf, job_tag, inv_type, inv_acc, n_grouped=q.get_num_node())
sol2 = inv2 * src
assert abs(sol.qnorm() - sol2.qnorm()) < 1e-20

# get_inv cached
inv3 = get_inv(gf, job_tag, inv_type, inv_acc, n_grouped=q.get_num_node())
inv4 = get_inv(gf, job_tag, inv_type, inv_acc, n_grouped=q.get_num_node())
assert inv3 is inv4

q.check_log_json(__file__)

q.timer_display()

qg.end_with_gpt()

q.displayln_info("CHECK: finished successfully.")
