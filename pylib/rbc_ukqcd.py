import qlat_gpt as qg
import qlat as q
import gpt as g
import gc
import rbc_ukqcd_params as rup

def get_total_site(job_tag : str):
    return rup.dict_params[job_tag]["total_site"]

@q.timer
def mk_sample_gauge_field(job_tag, fn):
    rs = q.RngState(f"seed {job_tag} {fn}").split("mk_sample_gauge_field")
    total_site = get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    gf = q.GaugeField(geo)
    gf.set_rand(rs, sigma = 0.25, n_step = 4)
    for i in range(4):
        q.gf_wilson_flow_step(gf, 0.05)
    gf.unitarize()
    return gf

@q.timer
def load_config(job_tag : str, fn : str):
    if not q.does_file_exist_sync_node(fn):
        if job_tag[:5] == "test-":
            return mk_sample_gauge_field(job_tag, fn)
        else:
            raise Exception(f"load_config '{fn}' does not exist.")
    gf = q.GaugeField()
    gf.load(fn)
    if (job_tag in rup.dict_params) and ("load_config_params" in rup.dict_params[job_tag]):
        params = rup.dict_params[job_tag]["load_config_params"]
        twist_boundary_at_boundary = params["twist_boundary_at_boundary"]
        for mu in range(4):
            lmom = twist_boundary_at_boundary[mu]
            if lmom != 0.0:
                q.displayln_info(f"load_config fn='{fn}' twist_boundary_at_boundary lmom={lmom} mu={mu}")
                gf.twist_boundary_at_boundary(lmom, mu)
    return gf

def get_fermion_params(job_tag, inv_type, inv_acc):
    return rup.dict_params[job_tag]["fermion_params"][inv_type][inv_acc]

def get_ls_from_fermion_params(fermion_params):
    if "omega" in fermion_params:
        return len(fermion_params["omega"])
    else:
        return fermion_params["Ls"]

def get_lanc_params(job_tag, inv_type, inv_acc = 0):
    return rup.dict_params[job_tag]["lanc_params"][inv_type][inv_acc]

def get_clanc_params(job_tag, inv_type, inv_acc = 0):
    return rup.dict_params[job_tag]["clanc_params"][inv_type][inv_acc]

@q.timer
def mk_eig(gf, job_tag, inv_type, inv_acc = 0):
    timer = q.Timer(f"py:mk_eig({job_tag},{inv_type},{inv_acc})", True)
    timer.start()
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    parity = g.odd
    params = get_lanc_params(job_tag, inv_type, inv_acc)
    q.displayln_info(f"mk_eig: job_tag={job_tag} inv_type={inv_type} inv_acc={inv_acc}")
    q.displayln_info(f"mk_eig: params={params}")
    fermion_params = params["fermion_params"]
    if "omega" in fermion_params:
        qm = g.qcd.fermion.zmobius(gpt_gf, fermion_params)
    else:
        qm = g.qcd.fermion.mobius(gpt_gf, fermion_params)
    w = g.qcd.fermion.preconditioner.eo2_ne(parity=parity)(qm)
    def make_src(rng):
        src = g.vspincolor(w.F_grid_eo)
        # src[:] = g.vspincolor([[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
        rng.cnormal(src)
        src.checkerboard(parity)
        return src
    pit = g.algorithms.eigen.power_iteration(**params["pit_params"])
    pit_ev, _, _ = pit(w.Mpc, make_src(g.random("lanc")));
    q.displayln_info(f"mk_eig: pit_ev={pit_ev}")
    #
    cheby = g.algorithms.polynomial.chebyshev(params["cheby_params"])
    irl = g.algorithms.eigen.irl(params["irl_params"])
    evec, ev = irl(cheby(w.Mpc), make_src(g.random("lanc")))
    evals = g.algorithms.eigen.evals(w.Mpc, evec, check_eps2=1e-6, real=True)
    g.mem_report()
    #
    timer.stop()
    return evec, evals

@q.timer
def mk_ceig(gf, job_tag, inv_type, inv_acc = 0):
    timer = q.Timer(f"py:mk_ceig({job_tag},{inv_type},{inv_acc})", True)
    timer.start()
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    parity = g.odd
    params = get_lanc_params(job_tag, inv_type, inv_acc)
    q.displayln_info(f"mk_ceig: job_tag={job_tag} inv_type={inv_type} inv_acc={inv_acc}")
    q.displayln_info(f"mk_ceig: params={params}")
    fermion_params = params["fermion_params"]
    if "omega" in fermion_params:
        qm = g.qcd.fermion.zmobius(gpt_gf, fermion_params)
    else:
        qm = g.qcd.fermion.mobius(gpt_gf, fermion_params)
    w = g.qcd.fermion.preconditioner.eo2_ne(parity=parity)(qm)
    def make_src(rng):
        src = g.vspincolor(w.F_grid_eo)
        # src[:] = g.vspincolor([[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
        rng.cnormal(src)
        src.checkerboard(parity)
        return src
    pit = g.algorithms.eigen.power_iteration(**params["pit_params"])
    pit_ev, _, _ = pit(w.Mpc, make_src(g.random("lanc")));
    q.displayln_info(f"mk_ceig: pit_ev={pit_ev}")
    #
    cheby = g.algorithms.polynomial.chebyshev(params["cheby_params"])
    irl = g.algorithms.eigen.irl(params["irl_params"])
    evec, ev = irl(cheby(w.Mpc), make_src(g.random("lanc")))
    evals = g.algorithms.eigen.evals(w.Mpc, evec, check_eps2=1e-6, real=True)
    g.mem_report()
    #
    inv = g.algorithms.inverter
    #
    cparams = get_clanc_params(job_tag, inv_type, inv_acc)
    q.displayln_info(f"mk_ceig: cparams={cparams}")
    #
    grid_coarse = g.block.grid(w.F_grid_eo, [ get_ls_from_fermion_params(fermion_params) ] + cparams["block"])
    nbasis = cparams["nbasis"]
    basis = evec[0:nbasis]
    b = g.block.map(grid_coarse, basis)
    for i in range(2):
        b.orthonormalize()
    del evec
    gc.collect()
    #
    ccheby = g.algorithms.polynomial.chebyshev(cparams["cheby_params"])
    cop = b.coarse_operator(ccheby(w.Mpc))
    #
    cstart = g.vcomplex(grid_coarse, nbasis)
    cstart[:] = g.vcomplex([1] * nbasis, nbasis)
    eps2 = g.norm2(cop * cstart - b.project * ccheby(w.Mpc) * b.promote * cstart) / g.norm2(cstart)
    g.message(f"Test coarse-grid promote/project cycle: {eps2}")
    cirl = g.algorithms.eigen.irl(cparams["irl_params"])
    cevec, cev = cirl(cop, cstart)
    #
    smoother = inv.cg(cparams["smoother_params"])(w.Mpc)
    smoothed_evals = []
    tmpf = g.lattice(basis[0])
    for i, cv in enumerate(cevec):
        tmpf @= smoother * b.promote * cv
        smoothed_evals = smoothed_evals + g.algorithms.eigen.evals(
            w.Mpc, [tmpf], check_eps2=10, real=True
        )
    g.mem_report()
    #
    timer.stop()
    return basis, cevec, smoothed_evals

@q.timer
def save_ceig(path, eig, job_tag, inv_type = 0, inv_acc = 0):
    if path is None:
        return
    save_params = get_clanc_params(job_tag, inv_type, inv_acc)["save_params"]
    nsingle = save_params["nsingle"]
    mpi = save_params["mpi"]
    fmt = g.format.cevec({"nsingle": nsingle, "mpi": [ 1 ] + mpi, "max_read_blocks": 8})
    g.save(path, eig, fmt);

@q.timer
def load_eig_lazy(path, job_tag, inv_type = 0, inv_acc = 0):
    # return ``None'' or a function ``load_eig''
    # ``load_eig()'' return the ``eig''
    if path is None:
        return None
    total_site = get_total_site(job_tag)
    fermion_params = get_fermion_params(job_tag, inv_type, inv_acc)
    grids = qg.get_fgrid(total_site, fermion_params)
    eig = None
    def load_eig():
        global eig
        if eig is None:
            eig = g.load(path, grids = grids)
        return eig
    return load_eig

@q.timer
def mk_gpt_inverter(gf, job_tag, inv_type, inv_acc, *,
        gt = None,
        mpi_split = None,
        n_grouped = 1,
        eig = None,
        eps = 1e-8,
        timer = True):
    if mpi_split is None:
        mpi_split = g.default.get_ivec("--mpi_split", None, 4)
        if mpi_split is not None:
            n_grouped = g.default.get_int("--grouped", 4)
    gpt_gf = qg.gpt_from_qlat(gf)
    pc = g.qcd.fermion.preconditioner
    if inv_type in [0, 1]:
        params = get_fermion_params(job_tag, inv_type, inv_acc)
        if eig is not None:
            # may need madwf
            params0 = get_fermion_params(job_tag, inv_type, inv_acc = 0)
            is_madwf = get_ls_from_fermion_params(params) != get_ls_from_fermion_params(params0)
        else:
            is_madwf = False
        if "omega" in params:
            qm = g.qcd.fermion.zmobius(gpt_gf, params)
        else:
            qm = g.qcd.fermion.mobius(gpt_gf, params)
        inv = g.algorithms.inverter
        if job_tag[:5] == "test-":
            cg_mp = inv.cg({"eps": eps, "maxiter": 100})
        elif inv_type == 0:
            cg_mp = inv.cg({"eps": eps, "maxiter": 200})
        elif inv_type == 1:
            cg_mp = inv.cg({"eps": eps, "maxiter": 300})
        else:
            raise Exception("mk_gpt_inverter")
        if mpi_split is None:
            cg_split = cg_mp
        else:
            cg_split = inv.split(cg_mp, mpi_split = mpi_split)
        if eig is not None:
            cg_defl = inv.coarse_deflate(eig[1], eig[0], eig[2])
            cg = inv.sequence(cg_defl, cg_split)
        else:
            cg = cg_split
        if inv_type == 0:
            slv_5d = inv.preconditioned(pc.eo2_ne(), cg)
        elif inv_type == 1:
            slv_5d = inv.preconditioned(pc.eo2_ne(), cg)
        else:
            raise Exception("mk_gpt_inverter")
        if is_madwf:
            gpt_gf_f = g.convert(gpt_gf, g.single)
            if "omega" in params0:
                qm0 = g.qcd.fermion.zmobius(gpt_gf_f, params0)
            else:
                qm0 = g.qcd.fermion.mobius(gpt_gf_f, params0)
            cg_pv_f = inv.cg({"eps": eps, "maxiter": 90})
            slv_5d_pv_f = inv.preconditioned(pc.eo2_ne(), cg_pv_f)
            slv_5d = pc.mixed_dwf(slv_5d, slv_5d_pv_f, qm0)
        if inv_acc == 0:
            maxiter = 1
        elif inv_acc == 1:
            maxiter = 2
        elif inv_acc == 2:
            maxiter = 200
        else:
            raise Exception("mk_gpt_inverter")
        slv_qm = qm.propagator(
                inv.defect_correcting(
                    inv.mixed_precision(
                        slv_5d, g.single, g.double),
                    eps=eps, maxiter=maxiter)).grouped(n_grouped)
        if timer is True:
            timer = q.Timer(f"py:inv({job_tag},{inv_type},{inv_acc})", True)
        elif timer is False:
            timer = q.TimerNone()
        inv_qm = qg.InverterGPT(inverter = slv_qm, timer = timer)
    else:
        raise Exception("mk_gpt_inverter")
    if gt is None:
        return inv_qm
    else:
        return q.InverterGaugeTransform(inverter = inv_qm, gt = gt)

@q.timer
def mk_qlat_inverter(gf, job_tag, inv_type, inv_acc, *, gt = None):
    timer = q.Timer(f"py:qinv({job_tag},{inv_type},{inv_acc})", True)
    if job_tag in ["24D", "32D"]:
        if inv_type == 0:
            fa = q.FermionAction(mass = 0.00107, m5 = 1.8, ls = 24, mobius_scale = 4.0)
            inv = q.InverterDomainWall(gf = gf, fa = fa, timer = timer)
            inv.set_stop_rsd(1e-8)
            inv.set_max_num_iter(200)
            if inv_acc == 0:
                maxiter = 1
            elif inv_acc == 1:
                maxiter = 2
            elif inv_acc == 2:
                maxiter = 50
            else:
                raise Exception("mk_qlat_inverter")
            inv.set_max_mixed_precision_cycle(maxiter)
        elif inv_type == 1:
            fa = q.FermionAction(mass = 0.0850, m5 = 1.8, ls = 24, mobius_scale = 4.0)
            inv = q.InverterDomainWall(gf = gf, fa = fa, timer = timer)
            inv.set_stop_rsd(1e-8)
            inv.set_max_num_iter(300)
            if inv_acc == 0:
                maxiter = 1
            elif inv_acc == 1:
                maxiter = 2
            elif inv_acc == 2:
                maxiter = 50
            else:
                raise Exception("mk_qlat_inverter")
            inv.set_max_mixed_precision_cycle(maxiter)
        else:
            raise Exception("mk_qlat_inverter")
    else:
        raise Exception("mk_qlat_inverter")
    if gt is None:
        return inv
    else:
        return q.InverterGaugeTransform(inverter = inv, gt = gt)

def mk_inverter(*args, **kwargs):
    return mk_gpt_inverter(*args, **kwargs)

@q.timer
def get_inv(gf, job_tag, inv_type, inv_acc, *, gt = None, mpi_split = None, n_grouped = 1, eig = None, eps = 1e-8, timer = True):
    tag = f"rbc_ukqcd.get_inv {id(gf)} {job_tag} {inv_type} {inv_acc} {id(gt)} {mpi_split} {n_grouped} {id(eig)} {eps} {id(timer)}"
    if tag in q.cache_inv:
        return q.cache_inv[tag]
    inv = mk_inverter(gf, job_tag, inv_type, inv_acc,
            gt = gt,
            mpi_split = mpi_split,
            n_grouped = n_grouped,
            eig = eig,
            eps = eps,
            timer = timer)
    q.cache_inv[tag] = inv
    return inv
