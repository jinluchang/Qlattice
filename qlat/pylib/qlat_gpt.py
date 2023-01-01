import qlat as q
import gpt as g
import math

import textwrap

def mk_grid(geo = None):
    if geo is None:
        l_size = 32 * 9 * 5
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
    return ctype.name + "," + str(list(total_site)) + "," + str(multiplicity) + "," + tag

def mk_gpt_field(ctype, geo):
    if ctype is q.ElemTypeColorMatrix:
        return g.mcolor(mk_grid(geo))
    elif ctype == q.ElemTypeWilsonMatrix:
        return g.mspincolor(mk_grid(geo))
    elif ctype == q.ElemTypeWilsonVector:
        return g.vspincolor(mk_grid(geo))
    elif ctype == q.ElemTypeComplex:
        return g.complex(mk_grid(geo))
    else:
        raise Exception("make_gpt_field")

@q.timer_verbose
def mk_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag):
    geo = q.Geometry(total_site, multiplicity)
    # q.displayln_info("qlat geo created")
    f_gpt = mk_gpt_field(ctype, geo)
    # q.displayln_info("gpt field made")
    f_qlat = q.Field(ctype, geo)
    # q.displayln_info("qlat field made")
    lexicographic_coordinates = g.coordinates(f_gpt)
    # q.displayln_info("gpt coordinates collected")
    buf = f_qlat.mview()
    # q.displayln_info("qlat mview made")
    if tag == "qlat_from_gpt":
        # q.displayln_info("qlat_from_gpt")
        qlat_from_gpt = g.copy_plan(buf, f_gpt)
        # q.displayln_info("plan initialized")
        qlat_from_gpt.destination += g.global_memory_view(
            f_gpt.grid,
            [[f_gpt.grid.processor, buf, 0, buf.nbytes]])
        # q.displayln_info("plan destination added")
        qlat_from_gpt.source += f_gpt.view[lexicographic_coordinates]
        # q.displayln_info("plan source added")
        qlat_from_gpt = qlat_from_gpt(local_only = True)
        # q.displayln_info("plan created")
        return qlat_from_gpt
    elif tag == "gpt_from_qlat":
        # q.displayln_info("gpt_from_qlat")
        gpt_from_qlat = g.copy_plan(f_gpt, buf)
        # q.displayln_info("plan initialized")
        gpt_from_qlat.source += g.global_memory_view(
            f_gpt.grid,
            [[f_gpt.grid.processor, buf, 0, buf.nbytes]])
        # q.displayln_info("plan source added")
        gpt_from_qlat.destination += f_gpt.view[lexicographic_coordinates]
        # q.displayln_info("plan destination added")
        gpt_from_qlat = gpt_from_qlat(local_only = True)
        # q.displayln_info("plan created")
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
    ctype = q.ElemTypeColorMatrix
    total_site = gpt_gf[0].grid.fdimensions
    multiplicity = 1
    tag = "qlat_from_gpt"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    geo = q.Geometry(total_site, 1)
    fs = [ q.FieldColorMatrix(geo) for i in range(4)]
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
    ctype = q.ElemTypeColorMatrix
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
    ctype = q.ElemTypeColorMatrix
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
    ctype = q.ElemTypeColorMatrix
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
    ctype = q.ElemTypeWilsonMatrix
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
    ctype = q.ElemTypeWilsonMatrix
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
    ctype = q.ElemTypeWilsonVector
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
    ctype = q.ElemTypeWilsonVector
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
    ctype = q.ElemTypeComplex
    total_site = gpt_fcs[0].grid.fdimensions
    multiplicity = 1
    tag = "qlat_from_gpt"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    geo = q.Geometry(total_site, 1)
    n = len(gpt_fcs)
    fs = [ q.FieldComplex(geo) for i in range(n) ]
    for i in range(n):
        plan(fs[i].mview(), gpt_fcs[i])
    if n == 1:
        return fs[0]
    ff = q.FieldComplex(q.Geometry(total_site, n))
    q.merge_fields(ff, fs)
    return ff

@q.timer
def gpt_from_qlat_complex(fc):
    assert isinstance(fc, q.FieldComplex)
    geo = fc.geo()
    ctype = q.ElemTypeComplex
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
    elif isinstance(obj, q.FieldComplex):
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
        assert isinstance(self.timer, (q.Timer, q.TimerNone,))
        assert isinstance(self.gpt_timer, (q.Timer, q.TimerNone,))

    def __mul__(self, prop_src):
        assert isinstance(prop_src, q.Prop) or isinstance(prop_src, q.FermionField4d) or isinstance(prop_src, list)
        self.timer.start()
        g_src = gpt_from_qlat(prop_src)
        g_sol = gpt_invert(g_src, self.inverter, self.gpt_timer)
        prop_sol = qlat_from_gpt(g_sol)
        self.timer.stop()
        return prop_sol

###

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

def line_search_quadratic(s, x, dx, dv0, df, step, *, max_c = 3):
    x = g.util.to_list(x)
    xp = g.copy(x)
    # ansatz: f(x) = a + b*(x-c)^2, then solve for c from dv1 and dv0
    # assume b > 0
    sv0 = g.group.inner_product(s, dv0)
    assert not math.isnan(sv0)
    sign = 1
    if sv0 == 0.0:
        return 0.0
    elif sv0 < 0:
        sign = -1
    c = 0.0
    sv_list = [ sv0, ]
    while True:
        dxp = []
        for dx_mu, s_mu in g.util.to_list(dx, s):
            mu = x.index(dx_mu)
            xp[mu] @= g(g.group.compose(sign * step * s_mu, xp[mu]))
            xp_mu = g.copy(xp[mu])
            g.project(xp[mu], "defect")
            project_diff2 = g.norm2(xp[mu] - xp_mu)
            if not (project_diff2 < 1e-8):
                g.message(f"line_search_quadratic: rank={g.rank()} project_diff={math.sqrt(project_diff2)} {sv_list}")
                if c == 0.0:
                    return None
                else:
                    return sign * c
            dxp.append(xp[mu])
        dv1 = df(xp, dxp)
        assert isinstance(dv1, list)
        sv1 = g.group.inner_product(s, dv1)
        sv_list.append(sv1)
        if math.isnan(sv1):
            g.message(f"line_search_quadratic: rank={g.rank()} {sv_list}")
            return None
        if sv0 > 0 and sv1 <= 0 or sv0 < 0 and sv1 >= 0:
            c += sv0 / (sv0 - sv1)
            return sign * c
        elif sv0 == 0.0:
            return sign * c
        else:
            c += 1
            sv0 = sv1
        if c > max_c:
            g.message(f"line_search_quadratic: rank={g.rank()} {sv_list}")
            return sign * c

class non_linear_cg(g.algorithms.base_iterative):

    @g.params_convention(
        eps = 1e-8,
        maxiter = 1000,
        step = 1e-3,
        log_functional_every = 10,
        line_search = line_search_quadratic,
        beta = g.algorithms.optimize.fletcher_reeves,
        max_c = 3,
    )
    def __init__(self, params):
        super().__init__()
        self.eps = params["eps"]
        self.maxiter = params["maxiter"]
        self.step = params["step"]
        self.nf = params["log_functional_every"]
        self.line_search = params["line_search"]
        self.beta = params["beta"]
        self.max_c = params["max_c"]

    def __call__(self, f):
        @self.timed_function
        def opt(x, dx, t):
            if self.maxiter <= 0:
                return False
            x = g.util.to_list(x)
            dx = g.util.to_list(dx)
            d_last = None
            s_last = None
            for i in range(self.maxiter):
                d = f.gradient(x, dx)
                assert isinstance(d, list)
                #
                if s_last is None:
                    beta = 0
                    s = d
                else:
                    beta = self.beta(d, d_last)
                    for nu in range(len(s)):
                        s[nu] = g(d[nu] + beta * s_last[nu])
                        if hasattr(s[nu].otype, "project"):
                            s[nu] = g.project(s[nu], "defect")
                #
                c = self.line_search(s, x, dx, d, f.gradient, -self.step, max_c = self.max_c)
                #
                rs = (
                    sum(g.norm2(d)) / sum([s.grid.gsites * s.otype.nfloats for s in d])
                ) ** 0.5
                #
                if c is None or math.isnan(c):
                    self.log(f"non_linear_cg: rank={g.rank()} c={c} reset s. iteration {i}: f(x) = {f(x):.15e}, |df|/sqrt(dof) = {rs:e}, beta = {beta}")
                    return False
                #
                for nu, x_mu in enumerate(dx):
                    x_mu @= g.group.compose(-self.step * c * s[nu], x_mu)
                    x_mu @= g.project(x_mu, "defect")
                #
                self.log_convergence(i, rs, self.eps)
                #
                if i % self.nf == 0:
                    self.log(
                        f"iteration {i}: f(x) = {f(x):.15e}, |df|/sqrt(dof) = {rs:e}, beta = {beta}, step = {c*self.step}"
                    )
                #
                if rs <= self.eps:
                    self.log(
                        f"converged in {i+1} iterations: f(x) = {f(x):.15e}, |df|/sqrt(dof) = {rs:e}"
                    )
                    return True
                #
                if abs(c) > self.max_c:
                    d_last = None
                    s_last = None
                    continue
                # keep last search direction
                d_last = d
                s_last = s
                #
            self.log(
                f"NOT converged in {i+1} iterations;  |df|/sqrt(dof) = {rs:e} / {self.eps:e}"
            )
            return False
            #
        return opt

###

@q.timer_verbose
def gauge_fix_coulomb(
        gf,
        *,
        gt = None,
        mpi_split = None,
        maxiter_gd = 10,
        maxiter_cg = 200,
        maxcycle_cg = 50000,
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
    #
    if gt is None:
        gtu = q.GaugeTransform(gf.geo())
        gtu.set_unit()
        V = gpt_from_qlat(gtu)
    else:
        V = gpt_from_qlat(gt)
    # load source
    U = gpt_from_qlat(gf)
    # split in time
    Nt = U[0].grid.gdimensions[3]
    g.message(f"Separate {Nt} time slices")
    Usep = [g.separate(u, 3) for u in U[0:3]]
    Vt = g.separate(V, 3)
    cache = {}
    if mpi_split is None:
        mpi_split = g.default.get_ivec("--mpi_split", [ 1, 1, 1, ], 3)
    split_grid = Usep[0][0].grid.split(mpi_split, Usep[0][0].grid.fdimensions)
    #
    g.message("Split grid")
    Usep_split = [g.split(Usep[mu], split_grid, cache) for mu in range(3)]
    Vt_split = g.split(Vt, split_grid, cache)
    #
    # optimizer
    opt = g.algorithms.optimize
    cg = non_linear_cg(
            maxiter = maxiter_cg,
            eps = eps,
            step = step,
            line_search = line_search_quadratic,
            log_functional_every = log_every,
            beta = opt.polak_ribiere,
            )
    gd = opt.gradient_descent(
            maxiter = maxiter_gd,
            eps = eps,
            step = step_gd,
            log_functional_every = log_every,
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
        if (rng is not None) and (gt is not None):
            rng.element(Vt_split[t])
        for i in range(maxcycle_cg):
            q.displayln(f"Running cg_cycle={i} local time slice {t} / {Nt_split} id_node={q.get_id_node()}")
            Vt_split[t] = g.project(Vt_split[t], "defect")
            if gd(fa)(Vt_split[t], Vt_split[t]):
                break
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
def check_gauge_fix_coulomb(gf, gt, eps = 1e-12):
    t_size = gf.geo().total_site()[3]
    V = gpt_from_qlat(gt)
    U = gpt_from_qlat(gf)
    Usep = [g.separate(u, 3) for u in U[0:3]]
    Vt = g.separate(V, 3)
    theta_list = []
    for t in range(t_size):
        f = g.qcd.gauge.fix.landau([Usep[mu][t] for mu in range(3)])
        dfv = f.gradient(Vt[t], Vt[t])
        theta = g.norm2(dfv).real / Vt[t].grid.gsites / dfv.otype.Nc
        theta_list.append(theta)
    s = sum(theta_list) / t_size
    q.displayln_info(f"check_gauge_fix_coulomb: theta slice average={s} eps={eps}")
    if not (s < eps):
        q.displayln_info(f"WARNING: check_gauge_fix_coulomb: failed with {theta_list}.")
    return s < eps
