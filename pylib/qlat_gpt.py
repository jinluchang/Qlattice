import gpt as g
import qlat as q

def mk_grid(geo = None):
    if geo is None:
        l_size = 8 * 3 * 5
        t_size = l_size * 2
        total_site = [l_size, l_size, l_size, t_size]
    else:
        total_site = geo.total_site()
    return g.grid(total_site, g.double)

def begin_with_gpt():
    grid = mk_grid()
    size_node = grid.mpi
    coor_node = grid.processor_coor
    id_node = q.index_from_coordinate(coor_node, size_node)
    q.begin(id_node, size_node)

def end_with_gpt():
    q.end()

def mk_qlat_gpt_copy_plan_key(ctype, total_site, multiplicity, tag):
    return ctype + "," + str(list(total_site)) + "," + str(multiplicity) + "," + tag

def mk_gpt_field(ctype, geo):
    if ctype == "ColorMatrix":
        return g.mcolor(mk_grid(geo))
    elif ctype == "WilsonMatrix":
        return g.mspincolor(mk_grid(geo))
    else:
        raise Exception("make_gpt_field")

@q.timer_verbose
def mk_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag):
    geo = q.Geometry(total_site, multiplicity)
    f_gpt = mk_gpt_field(ctype, geo)
    f_qlat = q.Field(ctype, geo)
    lexicographic_coordinates = g.coordinates(f_gpt)
    buf = f_qlat.mview()
    if tag == "qlat_from_gpt":
        qlat_from_gpt = g.copy_plan(buf, f_gpt)
        qlat_from_gpt.destination += g.global_memory_view(
            f_gpt.grid,
            [[f_gpt.grid.processor, buf, 0, buf.nbytes]])
        qlat_from_gpt.source += f_gpt.view[lexicographic_coordinates]
        qlat_from_gpt = qlat_from_gpt(local_only = True)
        return qlat_from_gpt
    elif tag == "gpt_from_qlat":
        gpt_from_qlat = g.copy_plan(f_gpt, buf)
        gpt_from_qlat.source += g.global_memory_view(
            f_gpt.grid,
            [[f_gpt.grid.processor, buf, 0, buf.nbytes]])
        gpt_from_qlat.destination += f_gpt.view[lexicographic_coordinates]
        gpt_from_qlat = gpt_from_qlat(local_only = True)
        return gpt_from_qlat
    else:
        q.displayln_info(tag)
        raise Exception("mk_qlat_gpt_copy_plan")

qlat_gpt_copy_plan_cache = {}

def get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag):
    key = mk_qlat_gpt_copy_plan_key(ctype, total_site, multiplicity, tag)
    if key in qlat_gpt_copy_plan_cache:
        return qlat_gpt_copy_plan_cache[key]
    else:
        plan = mk_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
        qlat_gpt_copy_plan_cache[key] = plan
        return plan

@q.timer
def qlat_from_gpt_gauge_field(gpt_gf):
    assert len(gpt_gf) == 4
    ctype = "ColorMatrix"
    total_site = gpt_gf[0].grid.fdimensions
    multiplicity = 1
    tag = "qlat_from_gpt"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    geo = q.Geometry(total_site, 1)
    fs = [ q.Field(ctype, geo) for i in range(4)]
    assert len(fs) == 4
    for i in range(4):
        plan(fs[i].mview(), gpt_gf[i])
    gf = q.GaugeField()
    q.merge_fields(gf, fs)
    return gf

@q.timer
def gpt_from_qlat_gauge_field(gf):
    assert isinstance(gf, q.GaugeField)
    geo = gf.geo()
    ctype = "ColorMatrix"
    total_site = geo.total_site()
    multiplicity = 1
    tag = "gpt_from_qlat"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    fs = [None] * 4
    q.split_fields(fs, gf)
    assert len(fs) == 4
    grid = mk_grid(geo)
    gpt_gf = [g.mcolor(grid) for i in range(4)]
    for i in range(4):
        plan(gpt_gf[i], fs[i].mview())
    return gpt_gf

@q.timer
def qlat_from_gpt_prop4d(gpt_prop):
    ctype = "WilsonMatrix"
    total_site = gpt_prop.grid.fdimensions
    multiplicity = 1
    tag = "qlat_from_gpt"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    geo = q.Geometry(total_site, 1)
    prop_msc = q.Prop(geo)
    plan(prop_msc.mview(), gpt_prop)
    prop_wm = q.Prop(geo)
    q.convert_wm_from_mspincolor_prop(prop_wm, prop_msc)
    return prop_wm

@q.timer
def gpt_from_qlat_prop4d(prop_wm):
    assert isinstance(prop_wm, q.Propagator4d)
    geo = prop_wm.geo()
    ctype = "WilsonMatrix"
    total_site = geo.total_site()
    multiplicity = 1
    tag = "gpt_from_qlat"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    prop_msc = q.Prop(geo)
    q.convert_mspincolor_from_wm_prop(prop_msc, prop_wm)
    grid = mk_grid(geo)
    gpt_prop = g.mspincolor(grid)
    plan(gpt_prop, prop_msc.mview())
    return gpt_prop

@q.timer
def qlat_from_gpt(gpt_obj):
    if repr(gpt_obj) == "lattice(ot_matrix_spin_color(4,3),double)":
        return qlat_from_gpt_prop4d(gpt_obj)
    elif repr(gpt_obj) == "[lattice(ot_matrix_su_n_fundamental_group(3),double), lattice(ot_matrix_su_n_fundamental_group(3),double), lattice(ot_matrix_su_n_fundamental_group(3),double), lattice(ot_matrix_su_n_fundamental_group(3),double)]":
        return qlat_from_gpt_gauge_field(gpt_obj)
    else:
        raise Exception("qlat_from_gpt")

@q.timer
def gpt_from_qlat(obj):
    if isinstance(obj, q.Propagator4d):
        return gpt_from_qlat_prop4d(obj)
    elif isinstance(obj, q.GaugeField):
        return gpt_from_qlat_gauge_field(obj)
    else:
        raise Exception("gpt_from_qlat")

@q.timer
def gpt_invert(src, inverter, timer = q.TimerNone()):
    sol = g.mspincolor(src.grid)
    timer.start()
    sol @= inverter * src
    timer.stop()
    return sol

class InverterGPT(q.Inverter):

    def __init__(self, *, inverter,
            timer = q.TimerNone(),
            gpt_timer = q.TimerNone()):
        self.inverter = inverter
        self.timer = timer
        self.gpt_timer = gpt_timer
        assert isinstance(self.timer, q.Timer)
        assert isinstance(self.gpt_timer, q.Timer)

    def __mul__(self, prop_src):
        assert isinstance(prop_src, q.Propagator4d)
        self.timer.start()
        g_src = gpt_from_qlat_prop4d(prop_src)
        g_sol = gpt_invert(g_src, self.inverter, self.gpt_timer)
        prop_sol = qlat_from_gpt_prop4d(g_sol)
        self.timer.stop()
        return prop_sol
