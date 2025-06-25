import qlat as q
import gc
import os
import pprint

from .rbc_ukqcd_params import get_param

def get_param_fermion(job_tag, inv_type, inv_acc):
    """
    Use param with lowewr `inv_acc` if the corresponding param does not exist.
    """
    while inv_acc >= 0:
        param = get_param(job_tag, "fermion_params", inv_type, inv_acc)
        if param is not None:
            return param
        inv_acc -= 1
    return None

def get_ls_from_fermion_params(fermion_params):
    if "omega" in fermion_params:
        return len(fermion_params["omega"])
    else:
        return fermion_params["Ls"]

def get_param_lanc(job_tag, inv_type, inv_acc=0):
    """
    Use param with lowewr `inv_acc` if the corresponding param does not exist.
    """
    while inv_acc >= 0:
        param = get_param(job_tag, "lanc_params", inv_type, inv_acc)
        if param is not None:
            return param
        inv_acc -= 1
    return None

def get_param_clanc(job_tag, inv_type, inv_acc=0):
    """
    Use param with lowewr `inv_acc` if the corresponding param does not exist.
    """
    while inv_acc >= 0:
        param = get_param(job_tag, "clanc_params", inv_type, inv_acc)
        if param is not None:
            return param
        inv_acc -= 1
    return None

def mk_pc_parity(job_tag, inv_type, inv_acc):
    import gpt as g
    return g.odd

def mk_pc_ne(job_tag, inv_type, inv_acc, *, eig=None, parity=None):
    import gpt as g
    fname = q.get_fname()
    params = get_param_fermion(job_tag, inv_type, inv_acc)
    if parity is None:
        parity = mk_pc_parity(job_tag, inv_type, inv_acc)
    pc = g.qcd.fermion.preconditioner
    if "omega" in params and eig is None:
        q.displayln_info(f"WARNING: {fname}: pc.eo2_kappa_ne() does not support split_cg. Try to avoid (use MDWF) if possible.")
        pc_ne = pc.eo2_kappa_ne(parity=parity)
    elif job_tag in [ "64I", "64I-pq", ]:
        pc_ne = pc.eo1_ne(parity=parity)
    else:
        pc_ne = pc.eo2_ne(parity=parity)
    return pc_ne

def mk_quark_matrix(job_tag, gpt_gf, inv_type, inv_acc):
    import gpt as g
    fermion_params = get_param_fermion(job_tag, inv_type, inv_acc)
    if "omega" in fermion_params:
        qm = g.qcd.fermion.zmobius(gpt_gf, fermion_params)
    else:
        qm = g.qcd.fermion.mobius(gpt_gf, fermion_params)
    return qm

@q.timer_verbose
def mk_eig(job_tag, gf, inv_type, inv_acc=0, *, parity=None, pc_ne=None):
    """
    Need get_param(job_tag, "lanc_params", inv_type, inv_acc)
    """
    import qlat_gpt as qg
    import gpt as g
    qtimer = q.Timer(f"py:mk_eig({job_tag},gf,{inv_type},{inv_acc})", True)
    qtimer.start()
    q.show_memory_usage()
    #
    if parity is None:
        parity = mk_pc_parity(job_tag, inv_type, inv_acc)
    if pc_ne is None:
        pc_ne = mk_pc_ne(job_tag, inv_type, inv_acc, parity=parity)
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    params = get_param_lanc(job_tag, inv_type, inv_acc)
    q.displayln_info(f"mk_eig: job_tag={job_tag} inv_type={inv_type} inv_acc={inv_acc}")
    q.displayln_info(f"mk_eig: params=\n{pprint.pformat(params)}")
    fermion_params = params["fermion_params"]
    assert fermion_params == get_param_fermion(job_tag, inv_type, inv_acc)
    qm = mk_quark_matrix(job_tag, gpt_gf, inv_type, inv_acc)
    w = pc_ne(qm)
    def make_src(rng):
        src = g.vspincolor(qm.F_grid_eo)
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
    evals, eps2 = g.algorithms.eigen.evals(w.Mpc, evec, real=True)
    #
    eig = qg.EigSystemGPT(evec=evec, evals=evals)
    q.show_memory_usage()
    qtimer.stop()
    return eig

@q.timer_verbose
def mk_ceig(job_tag, gf, inv_type, inv_acc=0, *, parity=None, pc_ne=None):
    """
    Need get_param(job_tag, "lanc_params", inv_type, inv_acc)
    Need get_param(job_tag, "clanc_params", inv_type, inv_acc)
    """
    import qlat_gpt as qg
    import gpt as g
    cparams = get_param_clanc(job_tag, inv_type, inv_acc)
    if cparams is None:
        return mk_eig(job_tag, gf, inv_type, inv_acc, parity=parity, pc_ne=pc_ne)
    qtimer = q.Timer(f"py:mk_ceig({job_tag},gf,{inv_type},{inv_acc})", True)
    qtimer.start()
    q.show_memory_usage()
    #
    if parity is None:
        parity = mk_pc_parity(job_tag, inv_type, inv_acc)
    if pc_ne is None:
        pc_ne = mk_pc_ne(job_tag, inv_type, inv_acc, parity=parity)
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    params = get_param_lanc(job_tag, inv_type, inv_acc)
    assert cparams["nbasis"] <= params["irl_params"]["Nstop"]
    fermion_params = params["fermion_params"]
    q.displayln_info(f"mk_ceig: job_tag={job_tag} inv_type={inv_type} inv_acc={inv_acc}")
    q.displayln_info(f"mk_ceig: params=\n{pprint.pformat(params)}")
    assert fermion_params == get_param_fermion(job_tag, inv_type, inv_acc)
    qm = mk_quark_matrix(job_tag, gpt_gf, inv_type, inv_acc)
    w = pc_ne(qm)
    def make_src(rng):
        src = g.vspincolor(qm.F_grid_eo)
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
    evals, eps2 = g.algorithms.eigen.evals(w.Mpc, evec, real=True)
    #
    inv = g.algorithms.inverter
    #
    q.displayln_info(f"mk_ceig: cparams=\n{pprint.pformat(cparams)}")
    #
    grid_coarse = g.block.grid(qm.F_grid_eo, [ get_ls_from_fermion_params(fermion_params) ] + cparams["block"])
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
        evals = g.algorithms.eigen.evals(
            w.Mpc, [tmpf], calculate_eps2=False, real=True
        )
        smoothed_evals = smoothed_evals + evals
    #
    eig = qg.EigSystemCompressedGPT(basis=basis, cevec=cevec, evals=smoothed_evals)
    q.show_memory_usage()
    qtimer.stop()
    return eig

@q.timer_verbose
def get_smoothed_evals(basis, cevec, gf, job_tag, inv_type, inv_acc=0, *, parity=None, pc_ne=None):
    """
    Need get_param(job_tag, "lanc_params", inv_type, inv_acc)
    Need get_param(job_tag, "clanc_params", inv_type, inv_acc)
    """
    import qlat_gpt as qg
    import gpt as g
    if parity is None:
        parity = mk_pc_parity(job_tag, inv_type, inv_acc)
    if pc_ne is None:
        pc_ne = mk_pc_ne(job_tag, inv_type, inv_acc, parity=parity)
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    params = get_param_lanc(job_tag, inv_type, inv_acc)
    cparams = get_param_clanc(job_tag, inv_type, inv_acc)
    assert cparams["nbasis"] <= params["irl_params"]["Nstop"]
    fermion_params = params["fermion_params"]
    assert fermion_params == get_param_fermion(job_tag, inv_type, inv_acc)
    qm = mk_quark_matrix(job_tag, gpt_gf, inv_type, inv_acc)
    w = pc_ne(qm)
    inv = g.algorithms.inverter
    grid_coarse = g.block.grid(qm.F_grid_eo, [ get_ls_from_fermion_params(fermion_params) ] + cparams["block"])
    b = g.block.map(grid_coarse, basis)
    smoother = inv.cg(cparams["smoother_params"])(w.Mpc)
    smoothed_evals = []
    tmpf = g.lattice(basis[0])
    for i, cv in enumerate(cevec):
        tmpf @= smoother * b.promote * cv
        evals = g.algorithms.eigen.evals(
            w.Mpc, [tmpf], calculate_eps2=False, real=True
        )
        smoothed_evals = smoothed_evals + evals
    return smoothed_evals

@q.timer_verbose
def save_eig(path, eig, job_tag, inv_type=0, inv_acc=0):
    """
    Need get_param(job_tag, "clanc_params", inv_type, inv_acc)
    """
    import qlat_gpt as qg
    if path is None:
        return
    assert isinstance(eig, q.EigSystem)
    if isinstance(eig, qg.EigSystemCompressedGPT):
        clanc_params = get_param_clanc(job_tag, inv_type, inv_acc)
        save_params = clanc_params["save_params"]
        nsingle = save_params["nsingle"]
        mpi = save_params["mpi"]
        eig.save(path, nsingle=nsingle, mpi=mpi)
    elif isinstance(eig, qg.EigSystemGPT):
        eig.save(path)
    else:
        assert False

@q.timer_verbose
def load_eig_lazy(path, job_tag, inv_type=0, inv_acc=0):
    """
    return ``None'' or a function ``load_eig''
    ``load_eig()'' return the ``eig''
    #
    Need get_param(job_tag, "total_site")
    Need get_param(job_tag, "fermion_params", inv_type, inv_acc)
    """
    import qlat_gpt as qg
    import gpt as g
    if path is None:
        q.displayln_info(f"load_eig_lazy: path is '{path}'")
        return None
    b1 = q.does_file_exist_qar_sync_node(os.path.join(path, "metadata.txt"))
    b2 = q.does_file_exist_qar_sync_node(os.path.join(path, "eigen-values.txt"))
    b3 = q.does_file_exist_qar_sync_node(os.path.join(path, "00/0000000000.compressed"))
    b4 = q.does_file_exist_qar_sync_node(os.path.join(path, "00.zip"))
    b5 = q.does_file_exist_qar_sync_node(os.path.join(path, "global"))
    b6 = q.does_file_exist_qar_sync_node(os.path.join(path, "index"))
    b7 = q.does_file_exist_qar_sync_node(os.path.join(path, "index.crc32"))
    if b1 and b2 and (b3 or b4):
        q.displayln_info(0, f"load_eig_lazy: '{path}' compressed eig available.")
        @q.lazy_call
        @q.timer_verbose
        def load_eig():
            eig = qg.EigSystemCompressedGPT()
            total_site = q.Coordinate(get_param(job_tag, "total_site"))
            fermion_params = get_param_fermion(job_tag, inv_type, inv_acc)
            eig.load(path, total_site=total_site, fermion_params=fermion_params)
            return eig
    elif b5 and b6 and b7:
        q.displayln_info(0, f"load_eig_lazy: '{path}' eig available.")
        @q.lazy_call
        @q.timer_verbose
        def load_eig():
            eig = qg.EigSystemGPT()
            eig.load(path)
            return eig
    else:
        q.displayln_info(f"load_eig_lazy: '{path}' {((b1, b2, b3, b4,), (b5, b6, b7,))}.")
        return None
    #
    return load_eig

def get_param_cg_mp_maxiter(job_tag, inv_type, inv_acc):
    maxiter = get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")
    if maxiter is not None:
        return maxiter
    if job_tag[:5] == "test-":
        maxiter = 50
    elif inv_acc >= 2:
        maxiter = 500
    elif inv_type == 0:
        maxiter = 200
    elif inv_type >= 1:
        maxiter = 300
    else:
        maxiter = 200
    return maxiter

@q.timer_verbose
def mk_gpt_inverter(
        gf, job_tag, inv_type, inv_acc,
        *,
        gt=None,
        mpi_split=None,
        n_grouped=None,
        eig=None,
        eps=1e-8,
        parity=None,
        pc_ne=None,
        qtimer=True,
        ):
    """
    Need get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")
    Need get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")
    Need get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "pv_maxiter", default=150) if is_madwf
    Need get_param(job_tag, "fermion_params", inv_type, inv_acc)
    Need get_param(job_tag, "fermion_params", inv_type, inv_acc=0) if eig is not None
    """
    import qlat_gpt as qg
    import gpt as g
    if mpi_split is None:
        mpi_split = g.default.get_ivec("--mpi_split", None, 4)
        if n_grouped is None:
            if mpi_split is not None:
                n_grouped = g.default.get_int("--grouped", 12)
            else:
                n_grouped = g.default.get_int("--grouped", 1)
    if n_grouped is None:
        n_grouped = 12
    q.displayln_info(f"mk_gpt_inverter: job_tag={job_tag} inv_type={inv_type} inv_acc={inv_acc} mpi_split={mpi_split} n_grouped={n_grouped}")
    gpt_gf = qg.gpt_from_qlat(gf)
    pc = g.qcd.fermion.preconditioner
    params = get_param_fermion(job_tag, inv_type, inv_acc)
    assert params is not None
    if eig is not None:
        # may need madwf
        params0 = get_param_fermion(job_tag, inv_type, inv_acc=0)
        is_madwf = get_ls_from_fermion_params(params) != get_ls_from_fermion_params(params0)
        if is_madwf and inv_type == 1:
            # do not use madwf for strange even when eig is available (do not use eig for exact solve)
            is_madwf = False
            eig = None
    else:
        is_madwf = False
    q.displayln_info(f"mk_gpt_inverter: set qm params={params}")
    assert params == get_param_fermion(job_tag, inv_type, inv_acc)
    qm = mk_quark_matrix(job_tag, gpt_gf, inv_type, inv_acc)
    inv = g.algorithms.inverter
    cg_mp = inv.cg({"eps": eps, "maxiter": get_param_cg_mp_maxiter(job_tag, inv_type, inv_acc)})
    cg_pv_f = inv.cg({"eps": eps, "maxiter": get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "pv_maxiter", default=150)})
    if mpi_split is None or mpi_split == False:
        cg_split = cg_mp
        mpi_split = None
    else:
        cg_split = inv.split(cg_mp, mpi_split=mpi_split)
        cg_pv_f = inv.split(cg_pv_f, mpi_split=mpi_split)
    q.displayln_info(f"mk_gpt_inverter: mpi_split={mpi_split} n_grouped={n_grouped}")
    if eig is not None:
        assert isinstance(eig, q.EigSystem)
        if isinstance(eig, qg.EigSystemCompressedGPT):
            cg_defl = inv.coarse_deflate(eig.cevec, eig.basis, eig.evals)
        elif isinstance(eig, qg.EigSystemGPT):
            cg_defl = inv.deflate(eig.evec, eig.evals)
        else:
            assert False
        cg = inv.sequence(cg_defl, cg_split)
    else:
        cg = cg_split
    if parity is None:
        parity = mk_pc_parity(job_tag, inv_type, inv_acc)
    if pc_ne is None:
        pc_ne = mk_pc_ne(job_tag, inv_type, inv_acc, eig=eig, parity=parity)
    slv_5d = inv.preconditioned(pc_ne, cg)
    q.displayln_info(f"mk_gpt_inverter: deal with is_madwf={is_madwf}")
    if is_madwf:
        gpt_gf_f = g.convert(gpt_gf, g.single)
        assert params0 == get_param_fermion(job_tag, inv_type, inv_acc=0)
        qm0 = mk_quark_matrix(job_tag, gpt_gf_f, inv_type, inv_acc=0)
        slv_5d_pv_f = inv.preconditioned(pc.eo2_ne(parity=parity), cg_pv_f)
        slv_5d = pc.mixed_dwf(slv_5d, slv_5d_pv_f, qm0)
    if inv_acc == 0:
        maxcycle = get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle", default=1)
    elif inv_acc == 1:
        maxcycle = get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle", default=2)
    elif inv_acc == 2:
        maxcycle = get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle", default=200)
    elif inv_acc == 3:
        maxcycle = get_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle", default=200)
    else:
        raise Exception("mk_gpt_inverter")
    q.sync_node()
    q.displayln_info(f"mk_gpt_inverter: eps={eps} maxcycle={maxcycle}")
    slv_qm = qm.propagator(
            inv.defect_correcting(
                inv.mixed_precision(
                    slv_5d, g.single, g.double),
                eps=eps, maxiter=maxcycle)).grouped(n_grouped)
    q.displayln_info(f"mk_gpt_inverter: make inv_qm")
    if qtimer is True:
        qtimer = q.Timer(f"py:inv({job_tag},{inv_type},{inv_acc})", True)
    elif qtimer is False:
        qtimer = q.TimerNone()
    inv_qm = qg.InverterGPT(inverter=slv_qm, qtimer=qtimer)
    if gt is None:
        return inv_qm
    else:
        return q.InverterGaugeTransform(inverter=inv_qm, gt=gt)

@q.timer_verbose
def mk_qlat_inverter(gf, job_tag, inv_type, inv_acc, *, gt=None):
    qtimer = q.Timer(f"py:qinv({job_tag},{inv_type},{inv_acc})", True)
    if job_tag in ["24D", "32D"]:
        if inv_type == 0:
            fa = q.FermionAction(mass=0.00107, m5=1.8, ls=24, mobius_scale=4.0)
            inv = q.InverterDomainWall(gf=gf, fa=fa, qtimer=qtimer)
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
            fa = q.FermionAction(mass=0.0850, m5=1.8, ls=24, mobius_scale=4.0)
            inv = q.InverterDomainWall(gf=gf, fa=fa, qtimer=qtimer)
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
        return q.InverterGaugeTransform(inverter=inv, gt=gt)

def mk_inverter(*args, **kwargs):
    return mk_gpt_inverter(*args, **kwargs)

@q.timer
def get_inv(
        gf, job_tag, inv_type, inv_acc, *,
        gt=None,
        mpi_split=None,
        n_grouped=None,
        eig=None,
        eps=1e-8,
        pc_ne=None,
        qtimer=True,
        ):
    tag = f"rbc_ukqcd.get_inv gf={id(gf)} {job_tag} inv_type={inv_type} inv_acc={inv_acc} gt={id(gt)} mpi_split={mpi_split} n_grouped={n_grouped} eig={id(eig)} pc_ne={id(pc_ne)} eps={eps} qtimer={id(qtimer)}"
    if tag in q.cache_inv:
        return q.cache_inv[tag]["inv"]
    inv = mk_inverter(
            gf, job_tag, inv_type, inv_acc,
            gt=gt,
            mpi_split=mpi_split,
            n_grouped=n_grouped,
            eig=eig,
            eps=eps,
            pc_ne=pc_ne,
            qtimer=qtimer,
            )
    q.cache_inv[tag] = {
            "inv": inv,
            "gf": gf,
            "gt": gt,
            "eig": eig,
            "pc_ne": pc_ne,
            "qtimer": qtimer,
            }
    return inv
