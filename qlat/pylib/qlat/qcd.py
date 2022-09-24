import qlat.cqlat as c

from qlat_utils import *

from qlat.field import *
from qlat.field_utils import *
from qlat.utils_io import *

class GaugeField(Field):

    def __init__(self, geo = None, *, ctype = None, multiplicity = None):
        Field.__init__(self, "ColorMatrix", geo, 4)

    def copy(self, is_copying_data = True):
        f = GaugeField()
        if is_copying_data:
            f @= self
        return f

    @timer
    def save(self, path):
        mk_file_dirs_info(path)
        return c.save_gauge_field(self, path)

    @timer
    def load(self, path):
        return c.load_gauge_field(self, path)

    def set_rand(self, rng, sigma = 0.5, n_step = 1):
        set_g_rand_color_matrix_field(self, rng, sigma, n_step)

    def unitarize(self):
        c.unitarize_color_matrix_field(self)

    def plaq(self):
        return gf_avg_plaq(self)

    def link_trace(self):
        return gf_avg_link_trace(self)

    def twist_boundary_at_boundary(self, lmom : float = -0.5, mu : int = 3):
        # modify in place
        gf_twist_boundary_at_boundary(self, lmom, mu)

    def show_info(self):
        gf_show_info(self)

###

class GaugeTransform(Field):

    def __init__(self, geo = None, *, ctype = None, multiplicity = None):
        Field.__init__(self, "ColorMatrix", geo, 1)

    def copy(self, is_copying_data = True):
        f = GaugeTransform()
        if is_copying_data:
            f @= self
        return f

    def set_rand(self, rng, sigma = 0.5, n_step = 1):
        set_g_rand_color_matrix_field(self, rng, sigma, n_step)

    def unitarize(self):
        c.unitarize_color_matrix_field(self)

    def __mul__(self, other):
        # other can be GaugeTransform, GaugeField, Prop, SelProp, PselProp, list
        from qlat.propagator import Prop, SelProp, PselProp
        if isinstance(other, GaugeTransform):
            gt = GaugeTransform()
            c.apply_gt_gt(gt, self, other)
            return gt
        elif isinstance(other, GaugeField):
            gf = GaugeField()
            c.apply_gt_gf(gf, self, other)
            return gf
        elif isinstance(other, Prop):
            prop = Prop()
            c.apply_gt_prop(prop, self, other)
            return prop
        elif isinstance(other, SelProp):
            prop = SelProp(other.fsel)
            c.apply_gt_sprop(prop, self, other)
            return prop
        elif isinstance(other, PselProp):
            prop = PselProp(other.psel)
            c.apply_gt_psprop(prop, self, other)
            return prop
        elif isinstance(other, list):
            return [ self * p for p in other ]
        else:
            raise Exception("GaugeTransform.__mul__")

    def inv(self):
        gt = GaugeTransform()
        c.gt_invert(gt, self)
        return gt

###

@timer_verbose
def gf_show_info(gf):
    assert isinstance(gf, GaugeField)
    displayln_info(f"gf_show_info: plaq = {gf.plaq():.16F} ; link_trace = {gf.link_trace():.16F}.")

def gf_avg_plaq(gf):
    assert isinstance(gf, GaugeField)
    return c.gf_avg_plaq(gf)

def gf_avg_link_trace(gf):
    assert isinstance(gf, GaugeField)
    return c.gf_avg_link_trace(gf)

def gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n = None):
    # wlf = Field("ColorMatrix", geo)
    # will only modify the m'th component
    # e.g. path = [ mu, mu, nu, -mu-1, -mu-1, ]
    # e.g. path = [ mu, nu, -mu-1, ], path_n = [ 2, 1, 2, ]
    if path_n is None:
        c.gf_wilson_line_no_comm(wlf, m, gf_ext, path)
    else:
        c.gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n)

def gf_wilson_lines_no_comm(gf_ext, path_list):
    # path_list = [ path_spec, ... ]
    # e.g. path_spec = [ mu, mu, nu, -mu-1, -mu-1, ]
    # e.g. path_spec = ([ mu, nu, -mu-1, ], [ 2, 1, 2, ],)
    # return wlf
    multiplicity = len(path_list)
    geo = geo_reform(gf_ext.geo(), multiplicity)
    wlf = Field("ColorMatrix", geo)
    for m, p in enumerate(path_list):
        if isinstance(p, tuple) and len(p) == 2:
            path, path_n = p
            gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n)
        else:
            path = p
            gf_wilson_line_no_comm(wlf, m, gf_ext, path)
    return wlf

def gf_avg_wilson_loop_normalized_tr(gf, l, t):
    assert isinstance(gf, GaugeField)
    assert isinstance(l, int)
    assert isinstance(t, int)
    return c.gf_avg_wilson_loop_normalized_tr(gf, l, t)

def set_g_rand_color_matrix_field(fc, rng, sigma, n_steps = 1):
    assert isinstance(fc, Field) and fc.ctype == "ColorMatrix"
    assert isinstance(rng, RngState)
    return c.set_g_rand_color_matrix_field(fc, rng, sigma, n_steps)

def gf_twist_boundary_at_boundary(gf : GaugeField, lmom : float = -0.5, mu : int = 3):
    # modify gf in place
    c.gf_twist_boundary_at_boundary(gf, lmom, mu)

def mk_left_expanded_gauge_field(gf):
    gf1 = field_expanded(gf, 1, 0)
    refresh_expanded_1(gf1)
    return gf1
