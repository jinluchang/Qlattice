import qlat_gpt as qg
import qlat as q
import gpt as g
import gc
import rbc_ukqcd_params as rup
import os
import pprint

from rbc_ukqcd_params import get_total_site
from rbc_ukqcd_params import mk_sample_gauge_field, load_config, load_config_lazy

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

@q.timer_verbose
def mk_eig(gf, job_tag, inv_type, inv_acc = 0):
    qtimer = q.Timer(f"py:mk_eig({job_tag},{inv_type},{inv_acc})", True)
    qtimer.start()
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    parity = g.odd
    params = get_lanc_params(job_tag, inv_type, inv_acc)
    q.displayln_info(f"mk_eig: job_tag={job_tag} inv_type={inv_type} inv_acc={inv_acc}")
    q.displayln_info(f"mk_eig: params=\n{pprint.pformat(params)}")
    fermion_params = params["fermion_params"]
    if "omega" in fermion_params:
        qm = g.qcd.fermion.zmobius(gpt_gf, fermion_params)
    else:
        qm = g.qcd.fermion.mobius(gpt_gf, fermion_params)
    w = g.qcd.fermion.preconditioner.eo2_ne(parity=parity)(qm)
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
    g.mem_report()
    #
    qtimer.stop()
    return evec, evals

@q.timer_verbose
def mk_ceig(gf, job_tag, inv_type, inv_acc = 0):
    qtimer = q.Timer(f"py:mk_ceig({job_tag},{inv_type},{inv_acc})", True)
    qtimer.start()
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    parity = g.odd
    params = get_lanc_params(job_tag, inv_type, inv_acc)
    cparams = get_clanc_params(job_tag, inv_type, inv_acc)
    assert cparams["nbasis"] <= params["irl_params"]["Nstop"]
    fermion_params = params["fermion_params"]
    q.displayln_info(f"mk_ceig: job_tag={job_tag} inv_type={inv_type} inv_acc={inv_acc}")
    q.displayln_info(f"mk_ceig: params=\n{pprint.pformat(params)}")
    if "omega" in fermion_params:
        qm = g.qcd.fermion.zmobius(gpt_gf, fermion_params)
    else:
        qm = g.qcd.fermion.mobius(gpt_gf, fermion_params)
    w = g.qcd.fermion.preconditioner.eo2_ne(parity=parity)(qm)
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
    g.mem_report()
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
            w.Mpc, [tmpf], calculate_eps2 = False, real=True
        )
        smoothed_evals = smoothed_evals + evals
    g.mem_report()
    #
    qtimer.stop()
    return basis, cevec, smoothed_evals

@q.timer_verbose
def get_smoothed_evals(basis, cevec, gf, job_tag, inv_type, inv_acc = 0):
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    parity = g.odd
    params = get_lanc_params(job_tag, inv_type, inv_acc)
    cparams = get_clanc_params(job_tag, inv_type, inv_acc)
    assert cparams["nbasis"] <= params["irl_params"]["Nstop"]
    fermion_params = params["fermion_params"]
    if "omega" in fermion_params:
        qm = g.qcd.fermion.zmobius(gpt_gf, fermion_params)
    else:
        qm = g.qcd.fermion.mobius(gpt_gf, fermion_params)
    w = g.qcd.fermion.preconditioner.eo2_ne(parity=parity)(qm)
    inv = g.algorithms.inverter
    grid_coarse = g.block.grid(qm.F_grid_eo, [ get_ls_from_fermion_params(fermion_params) ] + cparams["block"])
    b = g.block.map(grid_coarse, basis)
    smoother = inv.cg(cparams["smoother_params"])(w.Mpc)
    smoothed_evals = []
    tmpf = g.lattice(basis[0])
    for i, cv in enumerate(cevec):
        tmpf @= smoother * b.promote * cv
        evals = g.algorithms.eigen.evals(
            w.Mpc, [tmpf], calculate_eps2 = False, real=True
        )
        smoothed_evals = smoothed_evals + evals
    return smoothed_evals

@q.timer_verbose
def save_ceig(path, eig, job_tag, inv_type = 0, inv_acc = 0):
    if path is None:
        return
    save_params = get_clanc_params(job_tag, inv_type, inv_acc)["save_params"]
    nsingle = save_params["nsingle"]
    mpi = save_params["mpi"]
    fmt = g.format.cevec({"nsingle": nsingle, "mpi": [ 1 ] + mpi, "max_read_blocks": 8})
    q.mk_file_dirs_info(path)
    g.save(path, eig, fmt);

@q.timer_verbose
def load_eig_lazy(path, job_tag, inv_type = 0, inv_acc = 0):
    # return ``None'' or a function ``load_eig''
    # ``load_eig()'' return the ``eig''
    if path is None:
        q.displayln_info(f"load_eig_lazy: path is '{path}'")
        return None
    if not q.does_file_exist_sync_node(os.path.join(path, "metadata.txt")):
        q.displayln_info(f"load_eig_lazy: '{path}' has not metadata.")
        return None
    if not q.does_file_exist_sync_node(os.path.join(path, "eigen-values.txt")):
        q.displayln_info(f"load_eig_lazy: '{path}' has not eigen-values.")
        return None
    if not q.does_file_exist_sync_node(os.path.join(path, "00/0000000000.compressed")):
        q.displayln_info(f"load_eig_lazy: '{path}' has not data file '00/0000000000.compressed'.")
        return None
    #
    @q.timer_verbose
    def load_eig():
        total_site = get_total_site(job_tag)
        fermion_params = get_fermion_params(job_tag, inv_type, inv_acc)
        grids = qg.get_fgrid(total_site, fermion_params)
        eig = g.load(path, grids = grids)
        return eig
    #
    return q.lazy_call(load_eig)

def get_cg_mp_maxiter(job_tag, inv_type, inv_acc):
    cg_params = rup.dict_params[job_tag].get(f"cg_params-{inv_type}-{inv_acc}")
    if cg_params is not None:
        maxiter = cg_params.get("maxiter")
        if maxiter is not None:
            return maxiter
    if job_tag[:5] == "test-":
        maxiter = 50
    elif inv_acc == 2:
        maxiter = 500
    elif inv_type == 0:
        maxiter = 200
    elif inv_type in [ 1, 2, ]:
        maxiter = 300
    else:
        maxiter = 100
        raise Exception("get_cg_mp_maxiter")
    return maxiter

@q.timer_verbose
def mk_gpt_inverter(gf, job_tag, inv_type, inv_acc, *,
        gt = None,
        mpi_split = None,
        n_grouped = 1,
        eig = None,
        eps = 1e-8,
        qtimer = True):
    if mpi_split is None:
        mpi_split = g.default.get_ivec("--mpi_split", None, 4)
        if mpi_split is not None:
            n_grouped = g.default.get_int("--grouped", 4)
    q.displayln_info(f"mk_gpt_inverter: job_tag={job_tag} inv_type={inv_type} inv_acc={inv_acc} mpi_split={mpi_split} n_grouped={n_grouped}")
    gpt_gf = qg.gpt_from_qlat(gf)
    pc = g.qcd.fermion.preconditioner
    if inv_type in [ 0, 1, 2, ]:
        params = get_fermion_params(job_tag, inv_type, inv_acc)
        if eig is not None:
            # may need madwf
            params0 = get_fermion_params(job_tag, inv_type, inv_acc = 0)
            is_madwf = get_ls_from_fermion_params(params) != get_ls_from_fermion_params(params0)
            if is_madwf and inv_type == 1:
                # do not use madwf for strange even when eig is available (do not use eig for exact solve)
                is_madwf = False
                eig = None
        else:
            is_madwf = False
        q.displayln_info(f"mk_gpt_inverter: set qm params={params}")
        if "omega" in params:
            qm = g.qcd.fermion.zmobius(gpt_gf, params)
        else:
            qm = g.qcd.fermion.mobius(gpt_gf, params)
        inv = g.algorithms.inverter
        cg_mp = inv.cg({"eps": eps, "maxiter": get_cg_mp_maxiter(job_tag, inv_type, inv_acc)})
        cg_pv_f = inv.cg({"eps": eps, "maxiter": 150})
        if mpi_split is None or mpi_split == False:
            cg_split = cg_mp
            n_grouped = 1
            mpi_split = None
        else:
            cg_split = inv.split(cg_mp, mpi_split = mpi_split)
            cg_pv_f = inv.split(cg_pv_f, mpi_split = mpi_split)
        q.displayln_info(f"mk_gpt_inverter: mpi_split={mpi_split} n_grouped={n_grouped}")
        if eig is not None:
            cg_defl = inv.coarse_deflate(eig[1], eig[0], eig[2])
            cg = inv.sequence(cg_defl, cg_split)
        else:
            cg = cg_split
        if "omega" in params and eig is None:
            pc_ne = pc.eo2_kappa_ne()
        else:
            pc_ne = pc.eo2_ne()
        if inv_type == 0:
            slv_5d = inv.preconditioned(pc_ne, cg)
        elif inv_type in [ 1, 2, ]:
            slv_5d = inv.preconditioned(pc_ne, cg)
        else:
            raise Exception("mk_gpt_inverter")
        q.displayln_info(f"mk_gpt_inverter: deal with is_madwf={is_madwf}")
        if is_madwf:
            gpt_gf_f = g.convert(gpt_gf, g.single)
            if "omega" in params0:
                qm0 = g.qcd.fermion.zmobius(gpt_gf_f, params0)
            else:
                qm0 = g.qcd.fermion.mobius(gpt_gf_f, params0)
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
        q.sync_node()
        q.displayln_info(f"mk_gpt_inverter: eps={eps} max_cycle={maxiter}")
        slv_qm = qm.propagator(
                inv.defect_correcting(
                    inv.mixed_precision(
                        slv_5d, g.single, g.double),
                    eps=eps, maxiter=maxiter)).grouped(n_grouped)
        q.displayln_info(f"mk_gpt_inverter: make inv_qm")
        if qtimer is True:
            qtimer = q.Timer(f"py:inv({job_tag},{inv_type},{inv_acc})", True)
        elif qtimer is False:
            qtimer = q.TimerNone()
        inv_qm = qg.InverterGPT(inverter = slv_qm, qtimer = qtimer)
    else:
        raise Exception("mk_gpt_inverter")
    if gt is None:
        return inv_qm
    else:
        return q.InverterGaugeTransform(inverter = inv_qm, gt = gt)

@q.timer_verbose
def mk_qlat_inverter(gf, job_tag, inv_type, inv_acc, *, gt = None):
    qtimer = q.Timer(f"py:qinv({job_tag},{inv_type},{inv_acc})", True)
    if job_tag in ["24D", "32D"]:
        if inv_type == 0:
            fa = q.FermionAction(mass = 0.00107, m5 = 1.8, ls = 24, mobius_scale = 4.0)
            inv = q.InverterDomainWall(gf = gf, fa = fa, qtimer = qtimer)
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
            inv = q.InverterDomainWall(gf = gf, fa = fa, qtimer = qtimer)
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
def get_inv(gf, job_tag, inv_type, inv_acc, *, gt = None, mpi_split = None, n_grouped = 1, eig = None, eps = 1e-8, qtimer = True):
    tag = f"rbc_ukqcd.get_inv gf={id(gf)} {job_tag} inv_type={inv_type} inv_acc={inv_acc} gt={id(gt)} mpi_split={mpi_split} n_grouped={n_grouped} eig={id(eig)} eps={eps} qtimer={id(qtimer)}"
    if tag in q.cache_inv:
        return q.cache_inv[tag]["inv"]
    inv = mk_inverter(gf, job_tag, inv_type, inv_acc,
            gt = gt,
            mpi_split = mpi_split,
            n_grouped = n_grouped,
            eig = eig,
            eps = eps,
            qtimer = qtimer)
    q.cache_inv[tag] = {
            "inv": inv,
            "gf": gf,
            "gt": gt,
            "eig": eig,
            "qtimer": qtimer,
            }
    return inv
