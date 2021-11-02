from auto_contractor.eval import *
from auto_contractor.operators import *

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
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/vev-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer
def get_cexpr_meson_corr_psnk_psrc(vol):
    def calc_cexpr():
        exprs = [
                vol**2 * mk_pi_p("x2", True) * mk_pi_p("x1"),
                vol**2 * mk_k_p("x2", True) * mk_k_p("x1"),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        q.displayln_info(display_cexpr_raw(cexpr))
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_corr_psnk_psrc-cexpr.{vol}.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer
def get_cexpr_meson_corr_psnk_wsrc(vol):
    def calc_cexpr():
        exprs = [
                vol * mk_pi_p("x2", True) * mk_pi_p("t1"),
                vol * mk_k_p("x2", True) * mk_k_p("t1"),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        q.displayln_info(display_cexpr_raw(cexpr))
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_corr_psnk_wsrc-cexpr.{vol}.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

@q.timer
def get_cexpr_meson_corr_wsnk_wsrc():
    def calc_cexpr():
        exprs = [
                mk_pi_p("t2", True) * mk_pi_p("t1"),
                mk_k_p("t2", True) * mk_k_p("t1"),
                ]
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        q.displayln_info(display_cexpr_raw(cexpr))
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/meson_corr_wsnk_wsrc-cexpr.pickle")
    q.displayln_info(display_cexpr_raw(cexpr))
    return cexpr

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
        exprs_ops = exprs_odd_ops + exprs_even_ops
        exprs_src = [
            mk_k_0("t2_1") + "K0",
            mk_kpi_0_i1half("t2_1", "t2_2") + "Kpi_0_I1half",
            mk_kpi_0_i3halves("t2_1", "t2_2") + "Kpi_0_I3halves",
        ]
        exprs_snk = [
            mk_pipi_i20("t1_1", "t1_2", True) + "pipi_I2",
            mk_pipi_i0("t1_1", "t1_2", True) + "pipi_I0",
            mk_sigma("t1_1", True) + "sigma",
            mk_pi_0("t1_1", True) + "pi0",
            mk_expr(1) + "1",
        ]
        exprs_vac = [
            mk_pipi_i0("t1_1", "t1_2", True) + "pipi_I0",
            mk_sigma("t1_1", True) + "sigma",
        ]
        #    exprs_src_snk = [
        #        vol**2 * mk_
        #        ]
        exprs = []
        for expr_src in exprs_src:
            for expr_snk in exprs_snk:
                for expr_op in exprs_ops:
                    exprs.append(expr_snk * expr_op * expr_src)
        for expr_vac in exprs_vac:
            exprs.append(expr_vac)
        diagram_type_dict = dict()
        cexpr = contract_simplify_compile(*exprs, is_isospin_symmetric_limit = True, diagram_type_dict = diagram_type_dict)
        q.displayln_info(display_cexpr(cexpr))
        cexpr.collect_op()
        return cexpr
    cexpr = q.pickle_cache_call(calc_cexpr, f"cache/auto_contractor_cexpr/3f4f-cexpr.pickle")
    q.displayln_info(display_cexpr(cexpr))
    return cexpr

if __name__ == "__main__":
    get_cexpr_meson_corr_wsnk_wsrc()
    get_cexpr_vev()
    get_cexpr_3f4f_matching()
