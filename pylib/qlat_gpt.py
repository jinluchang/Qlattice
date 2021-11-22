import gpt as g
import qlat as q

import textwrap

def mk_grid(geo = None):
    if geo is None:
        l_size = 8 * 3 * 5
        t_size = l_size * 2
        total_site = [ l_size, l_size, l_size, t_size, ]
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
    elif ctype == "Complex":
        return g.complex(mk_grid(geo))
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
    fs = [ None, ] * 4
    q.split_fields(fs, gf)
    assert len(fs) == 4
    grid = mk_grid(geo)
    gpt_gf = [ None, ] * 4
    for i in range(4):
        gpt_gf[i] = g.mcolor(grid)
        plan(gpt_gf[i], fs[i].mview())
        fs[i] = None
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
    prop_wm = q.convert_wm_from_mspincolor(prop_msc)
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
    prop_msc = q.convert_mspincolor_from_wm(prop_wm)
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

@q.timer
def qlat_from_gpt_complex(gpt_fcs):
    assert isinstance(gpt_fcs, list)
    assert len(gpt_fcs) >= 1
    ctype = "Complex"
    total_site = gpt_fcs[0].grid.fdimensions
    multiplicity = 1
    tag = "qlat_from_gpt"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    geo = q.Geometry(total_site, 1)
    n = len(gpt_fcs)
    fs = [ q.Field(ctype, geo) for i in range(n) ]
    for i in range(n):
        plan(fs[i].mview(), gpt_fcs[i])
    if n == 1:
        return fs[0]
    ff = q.Field(ctype, q.Geometry(total_site, n))
    q.merge_fields(ff, fs)
    return ff

@q.timer
def gpt_from_qlat_complex(fc):
    assert isinstance(fc, q.Field)
    geo = fc.geo()
    ctype = "Complex"
    assert fc.ctype == ctype
    total_site = geo.total_site()
    multiplicity = 1
    tag = "gpt_from_qlat"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    n = geo.multiplicity()
    fs = [ None, ] * n
    q.split_fields(fs, fc)
    grid = mk_grid(geo)
    gpt_fcs = [ None, ] * n
    for i in range(n):
        gpt_fcs[i] = g.complex(grid)
        plan(gpt_fcs[i], fs[i].mview())
        fs[i] = None
    return gpt_fcs

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

def is_gpt_complex(obj):
    if isinstance(obj, list):
        for o in obj:
            if not (isinstance(o, g.core.lattice) and o.describe() == 'ot_complex_additive_group;none'):
                return False
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
    elif is_gpt_complex(gpt_obj):
        return qlat_from_gpt_complex(gpt_obj)
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
    elif isinstance(obj, q.Field) and obj.ctype == "Complex":
        return gpt_from_qlat_complex(obj)
    elif isinstance(obj, list):
        return [ gpt_from_qlat(p) for p in obj ]
    else:
        raise Exception("gpt_from_qlat")

@q.timer
def gpt_invert(src, inverter, qtimer = q.TimerNone()):
    qtimer.start()
    sol = g.eval(inverter * src)
    qtimer.stop()
    return sol

class InverterGPT(q.Inverter):

    def __init__(self, *, inverter,
            qtimer = q.TimerNone(),
            gpt_qtimer = q.TimerNone()):
        self.inverter = inverter
        self.timer = qtimer
        self.gpt_timer = gpt_qtimer
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

@q.timer_verbose
def save_gauge_field(gf, path):
    assert isinstance(gf, q.GaugeField)
    gpt_gf = gpt_from_qlat(gf)
    q.mk_file_dirs_info(path)
    g.save(path, gpt_gf, g.format.nersc())

@q.timer_verbose
def load_gauge_field(path):
    gpt_gf = g.load(path)
    return qlat_from_gpt(gpt_gf)

@q.timer_verbose
def gauge_fix_coulomb(
        gf,
        *,
        mpi_split = [ 1, 1, 1, ],
        maxiter_gd = 10,
        maxiter_cg = 200,
        maxcycle_cg = 50,
        log_every = 1,
        eps = 1e-12,
        step = 0.3,
        step_gd = 0.1,
        rng_seed = None,
        ):
    g.message(textwrap.dedent(f"""\
            Coulomb gauge fixer run with:
              mpi_split   = {mpi_split}
              maxiter_cg  = {maxiter_cg}
              maxiter_gd  = {maxiter_gd}
              maxcycle_cg = {maxcycle_cg}
              log_every   = {log_every}
              eps         = {eps}
              step        = {step}
              step_gd     = {step_gd}
              random      = {rng_seed}
            Note: convergence is only guaranteed for sufficiently small step parameter.
            """
            ))
    # create rng if needed
    rng = None if rng_seed is None else g.random(rng_seed)
    # load source
    U = gpt_from_qlat(gf)
    # split in time
    Nt = U[0].grid.gdimensions[3]
    g.message(f"Separate {Nt} time slices")
    Usep = [g.separate(u, 3) for u in U[0:3]]
    Vt = [g.mcolor(Usep[0][0].grid) for t in range(Nt)]
    cache = {}
    split_grid = Usep[0][0].grid.split(mpi_split, Usep[0][0].grid.fdimensions)
    #
    g.message("Split grid")
    Usep_split = [g.split(Usep[mu], split_grid, cache) for mu in range(3)]
    Vt_split = g.split(Vt, split_grid, cache)
    #
    # optimizer
    opt = g.algorithms.optimize
    cg = opt.non_linear_cg(
        maxiter=maxiter_cg,
        eps=eps,
        step=step,
        line_search=opt.line_search_quadratic,
        log_functional_every=log_every,
        beta=opt.polak_ribiere,
    )
    gd = opt.gradient_descent(
        maxiter=maxiter_gd,
        eps=eps,
        step=step_gd,
        log_functional_every=log_every,
    )
    #
    # Coulomb functional on each time-slice
    Nt_split = len(Vt_split)
    g.message(f"This rank has {Nt_split} time slices")
    #
    @q.timer
    def fix_t_slice(t):
        q.displayln(f"Run local time slice {t} / {Nt_split} id_node={q.get_id_node()}")
        f = g.qcd.gauge.fix.landau([Usep_split[mu][t] for mu in range(3)])
        fa = opt.fourier_accelerate.inverse_phat_square(Vt_split[t].grid, f)
        if rng is not None:
            rng.element(Vt_split[t])
        else:
            Vt_split[t] @= g.identity(Vt_split[t])
        if not gd(fa)(Vt_split[t], Vt_split[t]):
            for i in range(maxcycle_cg):
                Vt_split[t] = g.project(Vt_split[t], "defect")
                if cg(fa)(Vt_split[t], Vt_split[t]):
                    break
        q.displayln(f"Finish local time slice {t} / {Nt_split} id_node={q.get_id_node()}")
    #
    for t in range(Nt_split):
        fix_t_slice(t)
    #
    g.message("Unsplit")
    g.unsplit(Vt, Vt_split, cache)
    #
    g.message("Project to group (should only remove rounding errors)")
    Vt = [g.project(vt, "defect") for vt in Vt]
    #
    g.message("Test")
    # test results
    s = 0.0
    for t in range(Nt):
        f = g.qcd.gauge.fix.landau([Usep[mu][t] for mu in range(3)])
        dfv = f.gradient(Vt[t], Vt[t])
        theta = g.norm2(dfv).real / Vt[t].grid.gsites / dfv.otype.Nc
        s += theta
        g.message(f"theta[{t}] = {theta}")
        # g.message(f"V[{t}][0,0,0] = ", Vt[t][0, 0, 0])
    s /= Nt
    g.message(f"theta slice average = {s}")
    #
    # merge time slices
    V = g.merge(Vt, 3)
    gt = qlat_from_gpt(V)
    return gt

@q.timer
def check_gauge_fix_coulomb(gf, gt, eps = 1e-20):
    t_size = gf.geo().total_site()[3]
    V = gpt_from_qlat(gt)
    U = gpt_from_qlat(gf)
    Usep = [g.separate(u, 3) for u in U[0:3]]
    Vt = g.separate(V, 3)
    s = 0.0
    for t in range(t_size):
        f = g.qcd.gauge.fix.landau([Usep[mu][t] for mu in range(3)])
        dfv = f.gradient(Vt[t], Vt[t])
        theta = g.norm2(dfv).real / Vt[t].grid.gsites / dfv.otype.Nc
        s += theta
        q.displayln_info(f"theta[{t}] = {theta}")
    s /= t_size
    q.displayln_info(f"check_gauge_fix_coulomb: theta slice average = {s} eps={eps}")
    return s < eps
