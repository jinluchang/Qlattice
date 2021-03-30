import qlat_gpt as qg
import qlat as q
import gpt as g

def get_total_site(job_tag):
    if job_tag == "test-8nt16":
        return [8, 8, 8, 16]
    elif job_tag == "24D":
        return [24, 24, 24, 64]
    elif job_tag == "32D":
        return [32, 32, 32, 64]
    else:
        raise Exception("get_total_site")

def get_fermion_param(job_tag, inv_type, inv_accuracy):
    params = {
            "M5": 1.8,
            "b": 1.0,
            "c": 0.0,
            "boundary_phases": [1.0, 1.0, 1.0, 1.0],
            }
    if job_tag == "test-8nt16":
        params["b"] = 1.5
        params["c"] = 0.5
        if inv_type == 0:
            params["mass"] = 0.01
        elif inv_type == 1:
            params["mass"] = 0.04
        params["Ls"] = 8
    elif job_tag == "24D" or job_tag == "32D" or job_tag == "48D":
        params["b"] = 2.5
        params["c"] = 1.5
        if inv_type == 0:
            params["mass"] = 0.00107
        elif inv_type == 1:
            params["mass"] = 0.0850
        else:
            raise Exception("get_fermion_param")
        if inv_type == 0 and (inv_accuracy == 0 or inv_accuracy == 1):
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
        if inv_type == 0 and (inv_accuracy == 0 or inv_accuracy == 1):
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
        if inv_type == 0 and (inv_accuracy == 0 or inv_accuracy == 1):
            raise Exception("get_fermion_param")
        elif inv_type == 1 or inv_accuracy == 2:
            params["Ls"] = 24
        else:
            raise Exception("get_fermion_param")
    else:
        raise Exception("get_fermion_param")
    return params

@q.timer
def mk_gpt_inverter(gf, job_tag, inv_type, inv_accuracy, *, gt = None, n_grouped = 4):
    gpt_gf = qg.gpt_from_qlat(gf)
    pc = g.qcd.fermion.preconditioner
    if inv_type == 1:
        param = get_fermion_param(job_tag, inv_type, inv_accuracy)
        qm = None
        if "omega" in param:
            qm = g.qcd.fermion.zmobius(gpt_gf, param)
        else:
            qm = g.qcd.fermion.mobius(gpt_gf, param)
        inv = g.algorithms.inverter
        if inv_type == 0:
            cg_mp = inv.cg({"eps": 1e-8, "maxiter": 200})
        elif inv_type == 1:
            cg_mp = inv.cg({"eps": 1e-8, "maxiter": 300})
        else:
            raise Exception("mk_gpt_inverter")
        cg_split = inv.split(cg_mp, mpi_split = g.default.get_ivec("--mpi_split", None, 4))
        if inv_type == 0:
            slv_5d = inv.preconditioned(pc.eo2_ne(), cg_split)
        elif inv_type == 1:
            slv_5d = inv.preconditioned(pc.eo2_ne(), cg_split)
        else:
            raise Exception("mk_gpt_inverter")
        maxiter = 100
        if inv_accuracy == 0:
            maxiter = 1
        elif inv_accuracy == 1:
            maxiter = 2
        elif inv_accuracy == 2:
            maxiter = 50
        else:
            raise Exception("mk_gpt_inverter")
        slv_qm = qm.propagator(
                inv.defect_correcting(
                    inv.mixed_precision(
                        slv_5d, g.single, g.double),
                    eps=1e-8, maxiter=maxiter)).grouped(n_grouped)
        timer = q.Timer(f"py:inv({job_tag},{inv_type},{inv_accuracy})", True)
        inv_qm = qg.InverterGPT(inverter = slv_qm, timer = timer)
    elif inv_type == 0:
        raise Exception("mk_gpt_inverter")
    else:
        raise Exception("mk_gpt_inverter")
    if gt is None:
        return inv_qm
    else:
        return q.InverterGaugeTransform(inverter = inv_qm, gt = gt)

@q.timer
def mk_qlat_inverter(gf, job_tag, inv_type, inv_accuracy, *, gt = None):
    timer = q.Timer(f"py:qinv({job_tag},{inv_type},{inv_accuracy})", True)
    if job_tag == "24D" or job_tag == "32D":
        if inv_type == 0:
            fa = q.FermionAction(mass = 0.00107, m5 = 1.8, ls = 24, mobius_scale = 4.0)
            inv = q.InverterDomainWall(gf = gf, fa = fa, timer = timer)
            return inv
        elif inv_type == 1:
            fa = q.FermionAction(mass = 0.0850, m5 = 1.8, ls = 24, mobius_scale = 4.0)
            inv = q.InverterDomainWall(gf = gf, fa = fa, timer = timer)
            inv.set_stop_rsd(1e-8)
            inv.set_max_num_iter(300)
            maxiter = 100
            if inv_accuracy == 0:
                max_mixed_precision_cycle = -1
            elif inv_accuracy == 1:
                maxiter = 2
            elif inv_accuracy == 2:
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
