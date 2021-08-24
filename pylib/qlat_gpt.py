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
    elif ctype == "WilsonVector":
        return g.vspincolor(mk_grid(geo))
    else:
        raise Exception("make_gpt_field")

@q.timer_verbose
def mk_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag):
    geo = q.Geometry(total_site, multiplicity)
    q.displayln_info("qlat geo created")
    f_gpt = mk_gpt_field(ctype, geo)
    q.displayln_info("gpt field made")
    f_qlat = q.Field(ctype, geo)
    q.displayln_info("qlat field made")
    lexicographic_coordinates = g.coordinates(f_gpt)
    q.displayln_info("gpt coordinates collected")
    buf = f_qlat.mview()
    q.displayln_info("qlat mview made")
    if tag == "qlat_from_gpt":
        q.displayln_info("qlat_from_gpt")
        qlat_from_gpt = g.copy_plan(buf, f_gpt)
        q.displayln_info("plan initialized")
        qlat_from_gpt.destination += g.global_memory_view(
            f_gpt.grid,
            [[f_gpt.grid.processor, buf, 0, buf.nbytes]])
        q.displayln_info("plan destination added")
        qlat_from_gpt.source += f_gpt.view[lexicographic_coordinates]
        q.displayln_info("plan source added")
        qlat_from_gpt = qlat_from_gpt(local_only = True)
        q.displayln_info("plan created")
        return qlat_from_gpt
    elif tag == "gpt_from_qlat":
        q.displayln_info("gpt_from_qlat")
        gpt_from_qlat = g.copy_plan(f_gpt, buf)
        q.displayln_info("plan initialized")
        gpt_from_qlat.source += g.global_memory_view(
            f_gpt.grid,
            [[f_gpt.grid.processor, buf, 0, buf.nbytes]])
        q.displayln_info("plan source added")
        gpt_from_qlat.destination += f_gpt.view[lexicographic_coordinates]
        q.displayln_info("plan destination added")
        gpt_from_qlat = gpt_from_qlat(local_only = True)
        q.displayln_info("plan created")
        return gpt_from_qlat
    else:
        q.displayln_info(tag)
        raise Exception("mk_qlat_gpt_copy_plan")

cache_qlat_gpt_copy_plan = q.mk_cache("qlat_gpt_copy_plan")

def get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag):
    key = mk_qlat_gpt_copy_plan_key(ctype, total_site, multiplicity, tag)
    if key in cache_qlat_gpt_copy_plan:
        return cache_qlat_gpt_copy_plan[key]
    else:
        plan = mk_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
        cache_qlat_gpt_copy_plan[key] = plan
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
def qlat_from_gpt_gauge_transform(gpt_gt):
    ctype = "ColorMatrix"
    total_site = gpt_gt.grid.fdimensions
    multiplicity = 1
    tag = "qlat_from_gpt"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    geo = q.Geometry(total_site, 1)
    gt = q.GaugeTransform(geo)
    plan(gt.mview(), gpt_gt)
    return gt

@q.timer
def gpt_from_qlat_gauge_transform(gt):
    assert isinstance(gt, q.GaugeTransform)
    geo = gt.geo()
    ctype = "ColorMatrix"
    total_site = geo.total_site()
    multiplicity = 1
    tag = "gpt_from_qlat"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    grid = mk_grid(geo)
    gpt_gt = g.mcolor(grid)
    plan(gpt_gt, gt.mview())
    return gpt_gt

@q.timer
def qlat_from_gpt_prop(gpt_prop):
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
def gpt_from_qlat_prop(prop_wm):
    assert isinstance(prop_wm, q.Prop)
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
def qlat_from_gpt_ff4d(gpt_ff):
    ctype = "WilsonVector"
    total_site = gpt_ff.grid.fdimensions
    multiplicity = 1
    tag = "qlat_from_gpt"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    geo = q.Geometry(total_site, 1)
    ff = q.FermionField4d(geo)
    plan(ff.mview(), gpt_ff)
    return ff

@q.timer
def gpt_from_qlat_ff4d(ff):
    assert isinstance(ff, q.FermionField4d)
    geo = ff.geo()
    ctype = "WilsonVector"
    total_site = geo.total_site()
    multiplicity = 1
    tag = "gpt_from_qlat"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    grid = mk_grid(geo)
    gpt_ff = g.vspincolor(grid)
    plan(gpt_ff, ff.mview())
    return gpt_ff

def is_gpt_prop(obj):
    if isinstance(obj, g.core.lattice) and obj.describe() == "ot_matrix_spin_color(4,3);none":
        return True
    else:
        return False

def is_gpt_ff4d(obj):
    if isinstance(obj, g.core.lattice) and obj.describe() == "ot_vector_spin_color(4,3);none":
        return True
    else:
        return False

def is_gpt_gauge_field(obj):
    if isinstance(obj, list) and len(obj) == 4:
        for o in obj:
            if not (isinstance(o, g.core.lattice) and o.describe() == 'ot_matrix_su_n_fundamental_group(3);none'):
                return False
        return True
    else:
        return False

def is_gpt_gauge_transform(obj):
    if isinstance(obj, g.core.lattice) and obj.describe() == 'ot_matrix_su_n_fundamental_group(3);none':
        return True
    else:
        return False

@q.timer
def qlat_from_gpt(gpt_obj):
    if is_gpt_prop(gpt_obj):
        return qlat_from_gpt_prop(gpt_obj)
    elif is_gpt_gauge_transform(gpt_obj):
        return qlat_from_gpt_gauge_transform(gpt_obj)
    elif is_gpt_gauge_field(gpt_obj):
        return qlat_from_gpt_gauge_field(gpt_obj)
    elif is_gpt_ff4d(gpt_obj):
        return qlat_from_gpt_ff4d(gpt_obj)
    elif isinstance(gpt_obj, list):
        return [ qlat_from_gpt(p) for p in gpt_obj ]
    else:
        raise Exception("qlat_from_gpt")

@q.timer
def gpt_from_qlat(obj):
    if isinstance(obj, q.Prop):
        return gpt_from_qlat_prop(obj)
    elif isinstance(obj, q.GaugeTransform):
        return gpt_from_qlat_gauge_transform(obj)
    elif isinstance(obj, q.GaugeField):
        return gpt_from_qlat_gauge_field(obj)
    elif isinstance(obj, q.FermionField4d):
        return gpt_from_qlat_ff4d(obj)
    elif isinstance(obj, list):
        return [ gpt_from_qlat(p) for p in obj ]
    else:
        raise Exception("gpt_from_qlat")

@q.timer
def gpt_invert(src, inverter, timer = q.TimerNone()):
    timer.start()
    sol = g.eval(inverter * src)
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
        assert isinstance(prop_src, q.Prop) or isinstance(prop_src, q.FermionField4d) or isinstance(prop_src, list)
        self.timer.start()
        g_src = gpt_from_qlat(prop_src)
        g_sol = gpt_invert(g_src, self.inverter, self.gpt_timer)
        prop_sol = qlat_from_gpt(g_sol)
        self.timer.stop()
        return prop_sol

def get_fgrid(total_site, fermion_params):
    geo = q.Geometry(total_site, 1)
    gf = q.GaugeField(geo)
    gf.set_unit()
    gpt_gf = g.convert(gpt_from_qlat(gf), g.single)
    if "omega" in fermion_params:
        qm = g.qcd.fermion.zmobius(gpt_gf, fermion_params)
    else:
        qm = g.qcd.fermion.mobius(gpt_gf, fermion_params)
    return qm.F_grid_eo

def save_gauge_field(gf, path):
    assert isinstance(gf, q.GaugeField)
    gpt_gf = gpt_from_qlat(gf)
    g.save(path, gpt_gf, g.format.nersc())

def load_gauge_field(path):
    gpt_gf = g.load(path)
    return qlat_from_gpt(gpt_gf)
