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
