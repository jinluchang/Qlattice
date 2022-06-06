from auto_contractor.operators import *
from auto_contractor.eval import *
from jobs import *

@q.timer
def get_cexpr_vev():
    def calc_cexpr():
        s = new_spin_index()
        c = new_color_index()
        p = "x"
        exprs = [
                Qb("u", "x", s, c) * Qv("u", "x", s, c) + "u_bar*u",
                Qb("s", "x", s, c) * Qv("s", "x", s, c) + "s_bar*s",
                Qb("c", "x", s, c) * Qv("c", "x", s, c) + "c_bar*c",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contractor_cexpr/vev-cexpr")

@q.timer
def get_cexpr_meson_corr():
    def calc_cexpr():
        exprs = [
                mk_pi_p("x2", True) * mk_pi_p("x1") + "(pi   * pi)",
                mk_j5pi_mu("x2", 3) * mk_pi_p("x1") + "(a_pi * pi)",
                mk_k_p("x2", True)  * mk_k_p("x1")  + "(k    * k )",
                mk_j5k_mu("x2", 3)  * mk_k_p("x1")  + "(a_k  * k )",
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contractor_cexpr/meson_corr_cexpr")
    #cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_corr-cexpr.pickle")
    #q.displayln_info(display_cexpr_raw(cexpr))
    #return cexpr

@q.timer
def get_cexpr_meson_corr_with_env():
    def calc_cexpr():
        exprs_inner = [
                mk_pi_p("x2", True) * mk_pi_p("x1") + "(pi   * pi)",
                mk_j5pi_mu("x2", 3) * mk_pi_p("x1") + "(a_pi * pi)",
                mk_k_p("x2", True)  * mk_k_p("x1")  + "(k    * k )",
                mk_j5k_mu("x2", 3)  * mk_k_p("x1")  + "(a_k  * k )",
                ]
        exprs_outer = [
                mk_expr(1)                            + "vac_env ",
                mk_pi_p("x2p", True) * mk_pi_p("x1p") + "pi_p_env",
                mk_pi_0("x2p", True) * mk_pi_0("x1p") + "pi_0_env",
                mk_pi_m("x2p", True) * mk_pi_m("x1p") + "pi_m_env",
                ]
        exprs = []
        for e_o in exprs_outer:
            for e_i in exprs_inner:
                exprs.append(e_o * e_i)
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contractor_cexpr/meson_corr_with_env-cexpr")

@q.timer
def get_cexpr_corr_ww():
    def calc_cexpr():
        exprs = [
            mk_pi_p("t2", True) * mk_pi_p("t1") + "(pi   * pi)",
            mk_j5pi_mu("t2", 3) * mk_pi_p("t1") + "(a_pi * pi)",
            mk_k_p("t2", True)  * mk_k_p("t1")  + "(k    * k )",
            mk_j5k_mu("t2", 3)  * mk_k_p("t1")  + "(a_k  * k )",
            mk_k_0("t2", True)  * mk_k_0("t1")  + "(k  * k )",
            mk_kappa("t2", True)  * mk_kappa("t1")  + "(kappa  * kappa )",
            mk_Q0_b81("t2", "odd", True)  * mk_Q0_b81("t1", "odd")  + "(Q0o  * Q0o )",
            mk_Q1_b81("t2", "odd", True)  * mk_Q1_b81("t1", "odd")  + "(Q1o  * Q1o )",
            mk_Q2_b81("t2", "odd", True)  * mk_Q2_b81("t1", "odd")  + "(Q2o  * Q2o )",
            mk_Q3_b81("t2", "odd", True)  * mk_Q3_b81("t1", "odd")  + "(Q3o  * Q3o )",
            mk_Q4_b81("t2", "odd", True)  * mk_Q4_b81("t1", "odd")  + "(Q4o  * Q4o )",
            mk_Q5_b81("t2", "odd", True)  * mk_Q5_b81("t1", "odd")  + "(Q5o  * Q5o )",
            mk_Q6_b81("t2", "odd", True)  * mk_Q6_b81("t1", "odd")  + "(Q6o  * Q6o )",
            mk_Q7_b81("t2", "odd", True)  * mk_Q7_b81("t1", "odd")  + "(Q7o  * Q7o )",
            mk_Q8_b81("t2", "odd", True)  * mk_Q8_b81("t1", "odd")  + "(Q8o  * Q8o )",
            mk_Q0_b81("t2", "even", True)  * mk_Q0_b81("t1", "even")  + "(Q0e  * Q0e )",
            mk_Q1_b81("t2", "even", True)  * mk_Q1_b81("t1", "even")  + "(Q1e  * Q1e )",
            mk_Q2_b81("t2", "even", True)  * mk_Q2_b81("t1", "even")  + "(Q2e  * Q2e )",
            mk_Q3_b81("t2", "even", True)  * mk_Q3_b81("t1", "even")  + "(Q3e  * Q3e )",
            mk_Q4_b81("t2", "even", True)  * mk_Q4_b81("t1", "even")  + "(Q4e  * Q4e )",
            mk_Q5_b81("t2", "even", True)  * mk_Q5_b81("t1", "even")  + "(Q5e  * Q5e )",
            mk_Q6_b81("t2", "even", True)  * mk_Q6_b81("t1", "even")  + "(Q6e  * Q6e )",
            mk_Q7_b81("t2", "even", True)  * mk_Q7_b81("t1", "even")  + "(Q7e  * Q7e )",
            mk_Q8_b81("t2", "even", True)  * mk_Q8_b81("t1", "even")  + "(Q8e  * Q8e )",
        ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contractor_cexpr/meson_corr_cexpr")

@q.timer
def get_cexpr_3f4f_matching():
    @q.timer
    def calc_cexpr():
        exprs_odd_ops = [
            mk_Q0_b81("x", "odd") + "Q0_b81(o)",
            mk_Q1_b81("x", "odd") + "Q1_b81(o)",
            mk_Q2_b81("x", "odd") + "Q2_b81(o)",
            mk_Q3_b81("x", "odd") + "Q3_b81(o)",
            mk_Q4_b81("x", "odd") + "Q4_b81(o)",
            mk_Q5_b81("x", "odd") + "Q5_b81(o)",
            mk_Q6_b81("x", "odd") + "Q6_b81(o)",
            mk_Q7_b81("x", "odd") + "Q7_b81(o)",
            mk_Q8_b81("x", "odd") + "Q8_b81(o)",
        ]
        exprs_even_ops = [
            mk_Q0_b81("x", "even") + "Q0_b81(e)",
            mk_Q1_b81("x", "even") + "Q1_b81(e)",
            mk_Q2_b81("x", "even") + "Q2_b81(e)",
            mk_Q3_b81("x", "even") + "Q3_b81(e)",
            mk_Q4_b81("x", "even") + "Q4_b81(e)",
            mk_Q5_b81("x", "even") + "Q5_b81(e)",
            mk_Q6_b81("x", "even") + "Q6_b81(e)",
            mk_Q7_b81("x", "even") + "Q7_b81(e)",
            mk_Q8_b81("x", "even") + "Q8_b81(e)",
        ]
        #exprs_ops = exprs_odd_ops + exprs_even_ops
        exprs_src = [
            [ mk_k_0("t2_1") + "K0", -1, 0.5],
            [ mk_kpi_0_i1half("t2_1", "t2_2") + "Kpi_0_I1half_t2_2", 1, 0.5,],
            [ mk_kpi_0_i3halves("t2_1", "t2_2") + "Kpi_0_I3halves_t2_2", 1, 1.5,],
            [ mk_kpi_0_i1half("t2_1", "t2_3") + "Kpi_0_I1half_t2_3", 1, 0.5,],
            [ mk_kpi_0_i3halves("t2_1", "t2_3") + "Kpi_0_I3halves_t2_3", 1, 1.5,],
            [ mk_kpi_0_i1half("t2_1", "t2_4") + "Kpi_0_I1half_t2_4", 1, 0.5,],
            [ mk_kpi_0_i3halves("t2_1", "t2_4") + "Kpi_0_I3halves_t2_4", 1, 1.5,],
            [ mk_kappa("t2_1") + "kappa", 1, 0.5,],
        ]
        exprs_src_for2pt = [
            [ mk_Q0_b81("t2_1", "even", True) + "Q0_b81wsrc(e)", 1, 0.5,],
            [ mk_Q1_b81("t2_1", "even", True) + "Q1_b81wsrc(e)", 1, 0.5,],
            [ mk_Q2_b81("t2_1", "even", True) + "Q2_b81wsrc(e)", 1, 0.5,],
            [ mk_Q3_b81("t2_1", "even", True) + "Q3_b81wsrc(e)", 1, 0.5,],
            [ mk_Q4_b81("t2_1", "even", True) + "Q4_b81wsrc(e)", 1, 0.5,],
            [ mk_Q0_b81("t2_1", "odd", True) + "Q0_b81wsrc(o)", -1, 0.5,],
            [ mk_Q1_b81("t2_1", "odd", True) + "Q1_b81wsrc(o)", -1, 0.5,],
            [ mk_Q2_b81("t2_1", "odd", True) + "Q2_b81wsrc(o)", -1, 0.5,],
            [ mk_Q3_b81("t2_1", "odd", True) + "Q3_b81wsrc(o)", -1, 0.5,],
            [ mk_Q4_b81("t2_1", "odd", True) + "Q4_b81wsrc(o)", -1, 0.5,],
        ]
        exprs_snk = [
            [ mk_pipi_i20("t1_1", "t1_2", True) + "pipi_I2_t1_2", 1, 2,],
            [ mk_pipi_i0("t1_1", "t1_2", True) + "pipi_I0_t1_2", 1, 0,],
            [ mk_kk_i10("t1_1", "t1_2", True) + "KK_I1_t1_2", 1, 1,],
            [ mk_kk_i0("t1_1", "t1_2", True) + "KK_I0_t1_2", 1, 0,],
            [ mk_pipi_i20("t1_1", "t1_3", True) + "pipi_I2_t1_3", 1, 2,],
            [ mk_pipi_i0("t1_1", "t1_3", True) + "pipi_I0_t1_3", 1, 0,],
            [ mk_kk_i10("t1_1", "t1_3", True) + "KK_I1_t1_3", 1, 1,],
            [ mk_kk_i0("t1_1", "t1_3", True) + "KK_I0_t1_3", 1, 0,],
            [ mk_pipi_i20("t1_1", "t1_4", True) + "pipi_I2_t1_4", 1, 2,],
            [ mk_pipi_i0("t1_1", "t1_4", True) + "pipi_I0_t1_4", 1, 0,],
            [ mk_kk_i10("t1_1", "t1_4", True) + "KK_I1_t1_4", 1, 1,],
            [ mk_kk_i0("t1_1", "t1_4", True) + "KK_I0_t1_4", 1, 0,],
            [ mk_sigma("t1_1", True) + "sigma", 1, 0,],
            [ mk_pi_0("t1_1", True) + "pi0", -1, 1,],
            [ mk_a0_0("t1_1", True) + "a0_0", 1, 1,],
            [ mk_expr(1) + "1", 1, 0,],
        ]
        exprs_snk_for2pt = [
            [ mk_expr(1) + "1", 1, 0,],
        ]
        exprs_vac = [
            mk_pipi_i0("t1_1", "t1_2", True) + "pipi_I0_t1_2",
            mk_kk_i0("t1_1", "t1_2", True) + "KK_I0_t1_2",
            mk_pipi_i0("t1_1", "t1_3", True) + "pipi_I0_t1_3",
            mk_kk_i0("t1_1", "t1_3", True) + "KK_I0_t1_3",
            mk_pipi_i0("t1_1", "t1_4", True) + "pipi_I0_t1_4",
            mk_kk_i0("t1_1", "t1_4", True) + "KK_I0_t1_4",
            mk_sigma("t1_1", True) + "sigma",
        ]
        exprs_src_vec = [
            [
                mk_k_0_star_mu("t2_1",0) + "K0star0",
                mk_k_0_star_mu("t2_1",1) + "K0star1",
                mk_k_0_star_mu("t2_1",2) + "K0star2",
                -1, 0.5,
            ],
        ]
        exprs_snk_vec = [
            [
                mk_j10_mu("t1_1",0) + "rho0",
                mk_j10_mu("t1_1",1) + "rho1",
                mk_j10_mu("t1_1",2) + "rho2",
                -1, 1,
            ],
            [
                mk_j0_mu("t1_1",0) + "omega0",
                mk_j0_mu("t1_1",1) + "omega1",
                mk_j0_mu("t1_1",2) + "omega2",
                -1, 0,
            ],
        ]
        exprs = []
        for expr_src in exprs_src:
            for expr_snk in exprs_snk:
                if abs(expr_snk[2]-expr_src[2]) != 0.5 and abs(expr_snk[2]+expr_src[2]) != 0.5:
                    continue
                parity = expr_snk[1] * expr_src[1]
                exprs_ops = exprs_even_ops
                if parity == -1:
                    exprs_ops = exprs_odd_ops
                for expr_op in exprs_ops:
                    exprs.append(expr_snk[0] * expr_op * expr_src[0])
        for expr_src in exprs_src_for2pt:
            for expr_snk in exprs_snk_for2pt:
                if abs(expr_snk[2]-expr_src[2]) != 0.5 and abs(expr_snk[2]+expr_src[2]) != 0.5:
                    continue
                parity = expr_snk[1] * expr_src[1]
                exprs_ops = exprs_even_ops
                if parity == -1:
                    exprs_ops = exprs_odd_ops
                for expr_op in exprs_ops:
                    exprs.append(expr_snk[0] * expr_op * expr_src[0])
        for expr_src_vec in exprs_src_vec:
            for expr_snk_vec in exprs_snk_vec:
                if abs(expr_snk_vec[4]-expr_src_vec[4]) != 0.5 and abs(expr_snk_vec[4]+expr_src_vec[4]) != 0.5:
                    continue
                parity = expr_snk_vec[3] * expr_src_vec[3]
                exprs_ops = exprs_even_ops
                if parity == -1:
                    exprs_ops = exprs_odd_ops
                for expr_op in exprs_ops:
                    exprs.append(expr_snk_vec[0] * expr_op * expr_src_vec[0] + expr_snk_vec[1] * expr_op * expr_src_vec[1] + expr_snk_vec[2] * expr_op * expr_src_vec[2])
        for expr_vac in exprs_vac:
            exprs.append(expr_vac)
        diagram_type_dict = dict()
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    return cache_compiled_cexpr(calc_cexpr, f"cache/auto_contractor_cexpr/3f4f-cexpr")
    q.displayln_info(display_cexpr(cexpr))
    return cexpr

if __name__ == "__main__":
    get_cexpr_meson_corr()
    get_cexpr_meson_corr_with_env()
    get_cexpr_vev()
    get_cexpr_3f4f_matching()
