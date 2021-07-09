import qlat_gpt as qg
import qlat as q
import gpt as g
import gc

def get_total_site(job_tag):
    if job_tag == "test-4nt8":
        return [4, 4, 4, 8]
    elif job_tag == "test-4nt16":
        return [4, 4, 4, 16]
    elif job_tag == "test-8nt16":
        return [8, 8, 8, 16]
    elif job_tag == "test-16nt32":
        return [16, 16, 16, 32]
    elif job_tag == "test-32nt64":
        return [32, 32, 32, 64]
    elif job_tag == "test-48nt96":
        return [48, 48, 48, 96]
    elif job_tag == "test-64nt128":
        return [64, 64, 64, 128]
    elif job_tag == "test-96nt192":
        return [96, 96, 96, 192]
    elif job_tag == "test-128nt256":
        return [128, 128, 128, 256]
    elif job_tag == "16I":
        return [16, 16, 16, 32]
    elif job_tag in [ "24D", "24DH", "24I", ]:
        return [24, 24, 24, 64]
    elif job_tag == "32D":
        return [32, 32, 32, 64]
    elif job_tag == "48I":
        return [48, 48, 48, 96]
    else:
        raise Exception("get_total_site")

def get_fermion_param(job_tag, inv_type, inv_acc):
    params = {
            "M5": 1.8,
            "b": 1.0,
            "c": 0.0,
            "boundary_phases": [1.0, 1.0, 1.0, 1.0], # twist boundary after loading configuration
            }
    if job_tag[:5] == "test-":
        if inv_acc in [0, 1]:
            params["b"] = 1.5
            params["c"] = 0.5
            params["Ls"] = 8
        elif inv_acc in [2]:
            params["b"] = 10/8 + 0.5
            params["c"] = 10/8 - 0.5
            params["Ls"] = 10
        if inv_type == 0:
            params["mass"] = 0.01
        elif inv_type == 1:
            params["mass"] = 0.04
    elif job_tag in ["16I", "24I",]:
        params["b"] = 1.0
        params["c"] = 0.0
        params["Ls"] = 16
        if inv_type == 0:
            params["mass"] = 0.01
        elif inv_type == 1:
            params["mass"] = 0.04
    elif job_tag in ["24D", "32D", "48D"]:
        params["b"] = 2.5
        params["c"] = 1.5
        if inv_type == 0:
            params["mass"] = 0.00107
        elif inv_type == 1:
            params["mass"] = 0.0850
        else:
            raise Exception("get_fermion_param")
        if inv_type == 0 and (inv_acc == 0 or inv_acc == 1):
            params["omega"] = [
                    1.0903256131299373,
                    0.9570283702230611,
                    0.7048886040934104,
                    0.48979921782791747,
                    0.328608311201356,
                    0.21664245377015995,
                    0.14121112711957107,
                    0.0907785101745156,
                    0.05608303440064219 - 0.007537158177840385j,
                    0.05608303440064219 + 0.007537158177840385j,
                    0.0365221637144842 - 0.03343945161367745j,
                    0.0365221637144842 + 0.03343945161367745j
                    ]
        else:
            params["Ls"] = 24
    elif job_tag == "24DH":
        params["b"] = 2.5
        params["c"] = 1.5
        if inv_type == 0:
            params["mass"] = 0.0174
        elif inv_type == 1:
            params["mass"] = 0.0985
        else:
            raise Exception("get_fermion_param")
        if inv_type == 0 and (inv_acc == 0 or inv_acc == 1):
            params["omega"] = [
                    1.0903256131299373,
                    0.9570283702230611,
                    0.7048886040934104,
                    0.48979921782791747,
                    0.328608311201356,
                    0.21664245377015995,
                    0.14121112711957107,
                    0.0907785101745156,
                    0.05608303440064219 - 0.007537158177840385j,
                    0.05608303440064219 + 0.007537158177840385j,
                    0.0365221637144842 - 0.03343945161367745j,
                    0.0365221637144842 + 0.03343945161367745j
                    ]
        else:
            params["Ls"] = 24
    elif job_tag == "48I":
        params["b"] = 1.5
        params["c"] = 0.5
        if inv_type == 0:
            params["mass"] = 0.0006979
        elif inv_type == 1:
            params["mass"] = 0.03580
        else:
            raise Exception("get_fermion_param")
        if inv_type == 0 and (inv_acc == 0 or inv_acc == 1):
            raise Exception("get_fermion_param")
        elif inv_type == 1 or inv_acc == 2:
            params["Ls"] = 24
        else:
            raise Exception("get_fermion_param")
    else:
        raise Exception("get_fermion_param")
    return params

def get_ls_from_fermion_params(fermion_params):
    if "omega" in fermion_params:
        return len(fermion_params["omega"])
    else:
        return fermion_params["Ls"]

def get_lanc_params(job_tag, inv_type):
    inv_acc = 0
    fermion_params = get_fermion_param(job_tag, inv_type, inv_acc)
    pit_params = { "eps": 0.01, "maxiter": 500, "real": True }
    if job_tag[:5] == "test-":
        cheby_params = {"low": 0.05, "high": 5.5, "order": 50}
        irl_params = {
                "Nstop": 20,
                "Nk": 25,
                "Nm": 30,
                "resid": 1e-8,
                "betastp": 0.0,
                "maxiter": 20,
                "Nminres": 1,
                # "maxapply": 100
                }
    elif job_tag in ["16I"]:
        cheby_params = {"low": 0.005, "high": 5.5, "order": 100}
        irl_params = {
                "Nstop": 100,
                "Nk": 110,
                "Nm": 150,
                "resid": 1e-8,
                "betastp": 0.0,
                "maxiter": 20,
                "Nminres": 1,
                # "maxapply": 100
                }
    elif job_tag in ["24I"]:
        cheby_params = {"low": 0.001, "high": 5.5, "order": 100}
        irl_params = {
                "Nstop": 250,
                "Nk": 260,
                "Nm": 300,
                "resid": 1e-8,
                "betastp": 0.0,
                "maxiter": 20,
                "Nminres": 3,
                # "maxapply": 100
                }
    elif job_tag in ["24D"]:
        cheby_params = {"low": 2.97e-4, "high": 5.5, "order": 200}
        irl_params = {
                "Nstop": 1000,
                "Nk": 1040,
                "Nm": 1200,
                "resid": 1e-8,
                "betastp": 0.0,
                "maxiter": 20,
                "Nminres": 1,
                # "maxapply": 100
                }
    else:
        raise Exception("get_lanc_params")
    return {
            "fermion_params": fermion_params,
            "pit_params": pit_params,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            }

def get_clanc_params(job_tag, inv_type):
    inv_acc = 0
    if job_tag[:5] == "test-":
        block = [ 2, 2, 2, 2 ]
        nbasis = 20
        cheby_params = {"low": 0.20, "high": 5.5, "order": 50}
        irl_params = {
                "Nstop": 30,
                "Nk": 35,
                "Nm": 40,
                "resid": 1e-8,
                "betastp": 0.0,
                "maxiter": 20,
                "Nminres": 1,
                # "maxapply": 100
                }
        smoother_params = {"eps": 1e-6, "maxiter": 10}
    elif job_tag in ["16I"]:
        block = [ 2, 2, 2, 2 ]
        cheby_params = {"low": 0.015, "high": 5.5, "order": 100}
        nbasis = 100
        irl_params = {
                "Nstop": 300,
                "Nk": 310,
                "Nm": 350,
                "resid": 1e-8,
                "betastp": 0.0,
                "maxiter": 20,
                "Nminres": 1,
                # "maxapply": 100
                }
        smoother_params = {"eps": 1e-6, "maxiter": 10}
    elif job_tag in ["24I"]:
        block = [ 2, 2, 2, 2 ]
        cheby_params = {"low": 0.0025, "high": 5.5, "order": 100}
        nbasis = 250
        irl_params = {
                "Nstop": 500,
                "Nk": 510,
                "Nm": 550,
                "resid": 1e-8,
                "betastp": 0.0,
                "maxiter": 20,
                "Nminres": 4,
                #    "maxapply" : 100
                }
        smoother_params = {"eps": 1e-6, "maxiter": 10}
    elif job_tag in ["24D"]:
        block = [ 2, 3, 3, 2 ]
        cheby_params = {"low": 0.000684, "high": 5.5, "order": 200}
        nbasis = 1000
        irl_params = {
                "Nstop": 2000,
                "Nk": 2100,
                "Nm": 2600,
                "resid": 1e-8,
                "betastp": 0.0,
                "maxiter": 20,
                "Nminres": 1,
                # "maxapply": 100
                }
    else:
        raise Exception("get_clanc_params")
    return {
            "block": block,
            "nbasis": nbasis,
            "cheby_params": cheby_params,
            "irl_params": irl_params,
            "smoother_params": smoother_params,
            }

@q.timer
def mk_eig(gf, job_tag, inv_type):
    timer = q.Timer(f"py:mk_eig({job_tag},{inv_type})", True)
    timer.start()
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    parity = g.odd
    params = get_lanc_params(job_tag, inv_type)
    q.displayln_info(f"mk_eig: job_tag={job_tag} inv_type={inv_type}")
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
def mk_ceig(gf, job_tag, inv_type):
    timer = q.Timer(f"py:mk_ceig({job_tag},{inv_type})", True)
    timer.start()
    gpt_gf = g.convert(qg.gpt_from_qlat(gf), g.single)
    parity = g.odd
    params = get_lanc_params(job_tag, inv_type)
    q.displayln_info(f"mk_ceig: job_tag={job_tag} inv_type={inv_type}")
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
    cparams = get_clanc_params(job_tag, inv_type)
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
        param = get_fermion_param(job_tag, inv_type, inv_acc)
        if eig is not None:
            # may need madwf
            param0 = get_fermion_param(job_tag, inv_type, inv_acc = 0)
            is_madwf = get_ls_from_fermion_params(param) != get_ls_from_fermion_params(param0)
        else:
            is_madwf = False
        if "omega" in param:
            qm = g.qcd.fermion.zmobius(gpt_gf, param)
        else:
            qm = g.qcd.fermion.mobius(gpt_gf, param)
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
            if "omega" in param0:
                qm0 = g.qcd.fermion.zmobius(gpt_gf_f, param0)
            else:
                qm0 = g.qcd.fermion.mobius(gpt_gf_f, param0)
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
