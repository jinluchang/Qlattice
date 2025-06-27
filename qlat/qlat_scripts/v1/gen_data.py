from . import rbc_ukqcd as ru
from . import rbc_ukqcd_params as rup
from .jobs import *

# -----------------------------------------------------------------------------

@q.timer
def run_get_inverter(job_tag, traj, *, inv_type, get_gf, get_gt=None, get_eig=None):
    if None in [ get_gf, ]:
        return
    if get_gt is None:
        get_gt = lambda: None
    if get_eig is None:
        get_eig = lambda: None
    gf = get_gf()
    gt = get_gt()
    eig = get_eig()
    for inv_acc in [ 0, 1, 2, ]:
        ru.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_prop_wsrc_1(
        job_tag, traj,
        *,
        gf, gt, eig,
        idx, tslice,
        inv_type, inv_acc,
        ):
    q.check_stop()
    q.check_time_limit()
    fname = q.get_fname()
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    q.displayln_info(0, f"{fname}: {job_tag}/{traj} idx={idx} tag='{tag}'")
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    src = q.mk_wall_src(geo, tslice)
    sol = inv * src
    return sol

@q.timer
def compute_prop_wsrc_full(
        job_tag, traj,
        *,
        gf, gt, eig,
        idx, tslice,
        inv_type, inv_acc,
        sfw,
        ):
    fname = q.get_fname()
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in sfw:
        return None
    prop = compute_prop_wsrc_1(
            job_tag, traj,
            gf=gf, gt=gt, eig=eig,
            idx=idx, tslice=tslice,
            inv_type=inv_type, inv_acc=inv_acc,
            )
    q.qnorm_field(prop).save_double(sfw, tag + " ; qnorm_field")
    prop.save_double(sfw, tag, skip_if_exist=True)
    sfw.flush()

@q.timer_verbose
def compute_prop_wsrc_full_all(job_tag, traj, *,
                               inv_type, gf, gt, wi, eig):
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    path_s = f"{job_tag}/prop-wsrc-full-{inv_type_name}/traj-{traj}"
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    for inv_acc in [ 2, 1, ]:
        for p in wi:
            idx, tslice, inv_type_p, inv_acc_p=p
            if inv_type_p == inv_type and inv_acc_p == inv_acc:
                compute_prop_wsrc_full(
                        job_tag, traj,
                        gf=gf, gt=gt, eig=eig,
                        idx=idx, tslice=tslice, inv_type=inv_type, inv_acc=inv_acc,
                        sfw=sfw,
                        )
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

@q.timer(is_timer_fork=True)
def run_prop_wsrc_full(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt, get_wi):
    fname = q.get_fname()
    if None in [ get_gf, get_gt, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    path_f = f"{job_tag}/prop-wsrc-full-{inv_type_name}/traj-{traj}/geon-info.txt"
    path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}/geon-info.txt"
    if get_load_path(path_f) is not None:
        return
    if get_load_path(path_s) is not None:
        q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type_name} already have sparse wsrc prop. Skip calculating the full prop.")
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        wi = get_wi()
        compute_prop_wsrc_full_all(job_tag, traj,
                                   inv_type=inv_type, gf=gf, gt=gt, wi=wi,
                                   eig=eig)
        q.release_lock()
        return [ f"{fname} {job_tag} {traj} {inv_type} done", ]

# -----------------------------------------------------------------------------

@q.timer
def avg_weight_from_prop_full(geo, prop_nf_dict):
    """
    prop_nf_dict[(inv_type, tslice,)] => prop_nf, prop_nf_glb_sum_tslice
    type(prop_nf) = q.FieldRealD
    type(prop_nf_glb_sum_tslice) = q.SelectedPointsRealD
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    total_site = geo.total_site
    n_samples = [ 0, 0, ]
    avg_nf_glb_sum_tslice = [ 0, 0, ]
    for k, v in prop_nf_dict.items():
        inv_type, tslice = k
        prop_nf, prop_nf_glb_sum_tslice = v
        gst = np.roll(prop_nf_glb_sum_tslice[:].ravel(), -tslice)
        n_samples[inv_type] += 1
        avg_nf_glb_sum_tslice[inv_type] += gst
    for inv_type in [ 0, 1, ]:
        inv_type_name = inv_type_name_list[inv_type]
        assert n_samples[inv_type] == total_site[3]
        avg_nf_glb_sum_tslice[inv_type] = avg_nf_glb_sum_tslice[inv_type] / (geo.total_volume / total_site[3] * n_samples[inv_type])
        q.displayln_info(-1, fname, "avg_nf_glb_sum_tslice", inv_type_name, avg_nf_glb_sum_tslice[inv_type])
    local_tsize = geo.local_site[3]
    local_t_start = geo.coor_node[3] * local_tsize
    local_t_end = local_t_start + local_tsize
    local_tslices = slice(local_t_start, local_t_end, 1)
    local_field_shape = tuple(reversed(geo.local_site.to_list()))
    f_weight_avg = []
    for inv_type in [ 0, 1, ]:
        f = q.FieldRealD(geo, 1)
        q.set_zero(f)
        f_weight_avg.append(f)
    for k, v in prop_nf_dict.items():
        inv_type, tslice = k
        prop_nf, prop_nf_glb_sum_tslice = v
        avg_gst = avg_nf_glb_sum_tslice[inv_type][:].ravel()
        avg_gst_local = np.roll(avg_gst, tslice)[local_tslices]
        weight = prop_nf[:].reshape(local_field_shape) / avg_gst_local[:, None, None, None]
        # q.displayln_info(q.avg_err(weight.ravel()))
        view = f_weight_avg[inv_type][:].reshape(local_field_shape)
        view += weight / n_samples[inv_type]
    f_weight_final = q.FieldRealD(geo, 1)
    q.set_zero(f_weight_final)
    f_weight_final[:] = 0.5
    for inv_type in [ 0, 1, ]:
        f_weight_final[:] += f_weight_avg[inv_type][:] / 4
    return f_weight_avg, f_weight_final

@q.timer
def make_fsel_from_weight(f_weight, f_rand_01, rate):
    """
    f_weight is expected to be averaged around 1.
    """
    fname = q.get_fname()
    geo = f_weight.geo
    fsel = q.FieldSelection(geo)
    sel = f_weight[:].ravel() * rate >= f_rand_01[:].ravel()
    val = np.rint(f_weight[:].ravel()[sel] * 10**8).astype(int)
    assert np.all(val >= 0)
    fsel[sel] = val
    fsel.update()
    q.displayln_info(-1, f"{fname}: rate = {rate} ; expect_num = {geo.total_volume * rate} ; actual_num = {q.glb_sum(fsel.n_elems)}")
    return fsel

@q.timer
def make_psel_from_weight(f_weight, f_rand_01, rate):
    """
    f_weight is expected to be averaged around 1.
    """
    fsel = make_fsel_from_weight(f_weight, f_rand_01, rate)
    psel = fsel.to_psel()
    return psel

@q.timer_verbose
def compute_f_weight_from_wsrc_prop_full(job_tag, traj, *,
        fn_f_weight,
        inv_type_list,
        inv_type_name_list,
        path_f_list,
        fn_f_weight_type_list,
        ):
    fname = q.get_fname()
    assert get_load_path(fn_f_weight) is None
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    prop_nf_dict = dict()
    idx = 0
    inv_acc = 1
    for inv_type in inv_type_list:
        path_f = path_f_list[inv_type]
        sfr = q.open_fields(get_load_path(path_f), "r")
        available_tags = sfr.list()
        q.displayln_info(0, f"available_tags={available_tags}")
        for tslice in range(total_site[3]):
            tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
            if tag not in sfr:
                continue
            q.displayln_info(0, f"{fname}: idx={idx} tag='{tag}'")
            idx += 1
            prop_nf = q.FieldRealD(geo, 1)
            q.set_zero(prop_nf)
            prop_nf.load_double(sfr, tag + " ; qnorm_field")
            prop_nf_dict[(inv_type, tslice,)] = prop_nf, prop_nf.glb_sum_tslice()
    assert idx >= 1
    f_weight_avg, f_weight_final = avg_weight_from_prop_full(geo, prop_nf_dict)
    for inv_type in [ 0, 1, ]:
        inv_type_name = inv_type_name_list[inv_type]
        f_weight_avg[inv_type].save_double(get_save_path(fn_f_weight_type_list[inv_type]))
        q.displayln_info(-1, fname, "field-selection-weight", inv_type_name, f_weight_avg[inv_type].glb_sum_tslice()[:].ravel())
    f_weight_final.save_double(get_save_path(fn_f_weight))
    q.displayln_info(-1, fname, "field-selection-weight final", f_weight_final.glb_sum_tslice()[:].ravel())

@q.timer(is_timer_fork=True)
def run_f_weight_from_wsrc_prop_full(job_tag, traj):
    """
    return get_f_weight
        f_weight = get_f_weight()
    Or if wsrc_prop_full is not available
    return None
    #
    `f_weight` is of type `q.FieldRealD(geo, 1)`.
    #
    get_f_weight = run_f_weight_from_wsrc_prop_full(job_tag, traj)
    """
    fname = q.get_fname()
    fn_f_weight = f"{job_tag}/field-selection-weight/traj-{traj}/weight.field"
    @q.lazy_call
    @q.timer_verbose
    def get_f_weight():
        f_weight = q.FieldRealD()
        total_bytes = f_weight.load_double(get_load_path(fn_f_weight))
        assert total_bytes > 0
        return f_weight
    ret = get_f_weight
    if get_load_path(fn_f_weight) is not None:
        return ret
    inv_type_list = [ 0, 1, ]
    inv_type_name_list = [ "light", "strange", ]
    path_f_list = [
            f"{job_tag}/prop-wsrc-full-{inv_type_name}/traj-{traj}/geon-info.txt"
            for inv_type_name in inv_type_name_list
            ]
    fn_f_weight_type_list = [
            f"{job_tag}/field-selection-weight/traj-{traj}/weight-{inv_type_name}.field"
            for inv_type_name in inv_type_name_list
            ]
    for inv_type in inv_type_list:
        if get_load_path(path_f_list[inv_type]) is None:
            inv_type_name = inv_type_name_list[inv_type]
            q.displayln_info(-1, f"WARNING: {fname}: {job_tag} {traj} {inv_type_name} full prop is not available yet.")
            return None
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    compute_f_weight_from_wsrc_prop_full(
            job_tag, traj,
            fn_f_weight=fn_f_weight,
            inv_type_list=inv_type_list,
            inv_type_name_list=inv_type_name_list,
            path_f_list=path_f_list,
            fn_f_weight_type_list=fn_f_weight_type_list,
            )
    q.release_lock()
    return ret

@q.timer_verbose
def run_f_weight_uniform(job_tag, traj):
    """
    return get_f_weight
        f_weight = get_f_weight()
    #
    `f_weight` is of type `q.FieldRealD(geo, 1)`.
    #
    get_f_weight = run_f_weight_uniform(job_tag, traj)
    #
    Another option to obtain `get_f_weight` is using `run_f_weight_from_wsrc_prop_full(job_tag, traj)`.
    """
    fname = q.get_fname()
    fn_f_weight = f"{job_tag}/field-selection-weight/traj-{traj}/weight.field"
    @q.lazy_call
    @q.timer_verbose
    def get_f_weight():
        f_weight = q.FieldRealD()
        total_bytes = f_weight.load_double(get_load_path(fn_f_weight))
        assert total_bytes > 0
        return f_weight
    ret = get_f_weight
    if get_load_path(fn_f_weight) is not None:
        return ret
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    f_weight = q.FieldRealD(geo, 1)
    f_weight.set_unit()
    f_weight.save_double(get_save_path(fn_f_weight))
    q.release_lock()
    return ret

@q.timer_verbose
def run_f_weight_load(job_tag, traj):
    """
    return get_f_weight
        f_weight = get_f_weight()
    #
    `f_weight` is of type `q.FieldRealD(geo, 1)`.
    #
    get_f_weight = run_f_weight_load(job_tag, traj)
    #
    Just perform loading. Fail is data does not exist.
    """
    fname = q.get_fname()
    fn_f_weight = f"{job_tag}/field-selection-weight/traj-{traj}/weight.field"
    @q.lazy_call
    @q.timer_verbose
    def get_f_weight():
        f_weight = q.FieldRealD()
        total_bytes = f_weight.load_double(get_load_path(fn_f_weight))
        assert total_bytes > 0
        return f_weight
    ret = get_f_weight
    if get_load_path(fn_f_weight) is not None:
        return ret
    raise Exception(f"{fname}: 'fn_f_weight' does not exist.")

@q.timer_verbose
def run_f_rand_01(job_tag, traj):
    """
    return get_f_rand_01
        f_rand_01 = get_f_rand_01()
    #
    `f_rand_01` is of type `q.FieldRealD(geo, 1)`.
    #
    get_f_rand_01 = run_f_rand_01(job_tag, traj)
    """
    fname = q.get_fname()
    fn_f_rand_01 = f"{job_tag}/field-selection-weight/traj-{traj}/f-rand-01.field"
    @q.lazy_call
    @q.timer_verbose
    def get_f_rand_01():
        f_rand_01 = q.FieldRealD()
        total_bytes = f_rand_01.load_double(get_load_path(fn_f_rand_01))
        assert total_bytes > 0
        return f_rand_01
    ret = get_f_rand_01
    if get_load_path(fn_f_rand_01) is not None:
        return ret
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    seed = get_job_seed(job_tag)
    rs = q.RngState(f"{seed}-{traj}").split("run_sel_from_wsrc_prop_full").split("f_rand_01")
    f_rand_01 = q.FieldRealD(geo, 1)
    f_rand_01.set_rand(rs, 1.0, 0.0)
    f_rand_01.save_double(get_save_path(fn_f_rand_01))
    q.release_lock()
    return ret

@q.timer_verbose
def run_fsel_prob(job_tag, traj, *, get_f_rand_01, get_f_weight):
    """
    return get_fsel_prob
    #
        fsel_prob = get_fsel_prob()
        fsel = fsel_prob.fsel
    #
    One can set `get_f_rand_01=None` and `get_f_weight=None` if one intend to load existing data on `fsel_prob` instead of generating new `fsel_prob`.
    #
    `get_f_weight` is of type `lambda : q.FieldRealD(geo, 1)`.
    `get_f_rand_01` is of type `lambda : q.FieldRealD(geo, 1)`.
    #
    get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = lambda : get_fsel_prob().fsel
    #
    For old data format, in which case the `fsel-prob.sfield` data is absent, we assume uniform probability AND will combine the `fsel` with `psel` obtained from `run_psel_prob` function.
    """
    fname = q.get_fname()
    fn_fsel = f"{job_tag}/field-selection/traj-{traj}.field"
    fn_fsel_prob = f"{job_tag}/field-selection-weight/traj-{traj}/fsel-prob.sfield"
    @q.lazy_call
    @q.timer_verbose
    def get_fsel_prob():
        fsel = q.FieldSelection()
        total_size = fsel.load(get_load_path(fn_fsel))
        assert total_size > 0
        fsel_prob = q.SelectedFieldRealD(fsel, 1)
        total_size = fsel_prob.load_double(get_load_path(fn_fsel_prob))
        assert total_size > 0
        return fsel_prob
    ret = get_fsel_prob
    if get_load_path(fn_fsel) is not None:
        if get_load_path(fn_fsel_prob) is not None:
            return ret
        else:
            q.displayln_info(f"{fname}: WARNING: field-selection exist but prob is not available. Assuming loading load old data format.")
            get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
            @q.lazy_call
            @q.timer_verbose
            def get_fsel_prob_old():
                fsel = q.FieldSelection()
                total_size = fsel.load(get_load_path(fn_fsel))
                assert total_size > 0
                psel = get_psel_prob().psel
                fsel.add_psel(psel)
                total_volume = fsel.total_site.volume()
                fsel_prob = q.SelectedFieldRealD(fsel, 1)
                fsel_prob[:] = q.glb_sum(len(fsel)) / total_volume
                return fsel_prob
            return get_fsel_prob_old
    if get_f_rand_01 is None:
        q.displayln_info(-1, f"{fname}: get_f_rand_01 is None")
        return None
    if get_f_weight is None:
        q.displayln_info(-1, f"{fname}: get_f_weight is None")
        return None
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    fsel_rate = get_param(job_tag, "field_selection_fsel_rate")
    q.displayln_info(-1, fname, f"fsel_rate = {fsel_rate}")
    assert fsel_rate is not None
    assert get_load_path(fn_fsel) is None
    assert get_load_path(fn_fsel_prob) is None
    f_rand_01 = get_f_rand_01()
    f_weight = get_f_weight()
    fsel = make_fsel_from_weight(f_weight, f_rand_01, fsel_rate)
    fsel.save(get_save_path(fn_fsel))
    fsel_prob = q.SelectedFieldRealD(fsel, 1)
    fsel_prob @= f_weight
    fsel_prob[:] = np.minimum(1.0, fsel_prob[:] * fsel_rate)
    fsel_prob.save_double(get_save_path(fn_fsel_prob))
    q.release_lock()
    return ret

@q.timer_verbose
def run_psel_prob(job_tag, traj, *, get_f_rand_01, get_f_weight, tag=None):
    """
    return get_psel_prob
    #
        psel_prob = get_psel_prob()
        psel = psel_prob.psel
    #
    One can set `get_f_rand_01=None` and `get_f_weight=None` if one intend to load existing data on `psel_prob` instead of generating new `psel_prob`.
    #
    `get_f_weight` is of type `lambda : q.FieldRealD(geo, 1)`.
    `get_f_rand_01` is of type `lambda : q.FieldRealD(geo, 1)`.
    #
    tag can be "small", "median", "large", etc
    #
    get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel = lambda : get_psel_prob().psel
    """
    fname = q.get_fname()
    if tag is None:
        fn_psel = f"{job_tag}/points-selection/traj-{traj}.lati"
        fn_psel_prob = f"{job_tag}/field-selection-weight/traj-{traj}/psel-prob.lat"
    else:
        fn_psel = f"{job_tag}/psel_{tag}/traj-{traj}/psel.lati"
        fn_psel_prob = f"{job_tag}/psel_{tag}/traj-{traj}/psel-prob.lat"
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    @q.lazy_call
    @q.timer_verbose
    def get_psel_prob():
        psel = q.PointsSelection()
        psel.load(get_load_path(fn_psel))
        psel_prob = q.SelectedPointsRealD(psel, 1)
        psel_prob.load(get_load_path(fn_psel_prob))
        return psel_prob
    ret = get_psel_prob
    if get_load_path(fn_psel) is not None:
        if get_load_path(fn_psel_prob) is not None:
            return ret
        else:
            q.displayln_info(f"{fname}: WARNING: points-selection exist but prob is not available. Assuming loading load old data format.")
            total_volume = total_site.volume()
            @q.lazy_call
            @q.timer_verbose
            def get_psel_prob_old():
                psel = q.PointsSelection()
                psel.load(get_load_path(fn_psel))
                psel_prob = q.SelectedPointsRealD(psel, 1)
                psel_prob[:] = len(psel) / total_volume
                return psel_prob
            return get_psel_prob_old
    if get_f_rand_01 is None:
        q.displayln_info(-1, f"{fname}: get_f_rand_01 is None")
        return None
    if get_f_weight is None:
        q.displayln_info(-1, f"{fname}: get_f_weight is None")
        return None
    if tag is None:
        if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
            return None
    else:
        if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{tag}"):
            return None
    if tag is None:
        psel_rate = get_param(job_tag, "field_selection_psel_rate")
    else:
        psel_rate = get_param(job_tag, f"field_selection_psel_rate_{tag}")
    q.displayln_info(-1, fname, f"tag='{tag}' ; psel_rate={psel_rate} .")
    assert psel_rate is not None
    assert get_load_path(fn_psel) is None
    assert get_load_path(fn_psel_prob) is None
    f_rand_01 = get_f_rand_01()
    f_weight = get_f_weight()
    psel = make_psel_from_weight(f_weight, f_rand_01, psel_rate)
    psel.save(get_save_path(fn_psel))
    psel_prob = q.SelectedPointsRealD(psel, 1)
    psel_prob @= f_weight
    psel_prob[:] = np.minimum(1.0, psel_prob[:] * psel_rate)
    psel_prob.save(get_save_path(fn_psel_prob))
    q.release_lock()
    return ret

@q.timer
def run_fsel_from_fsel_prob(get_fsel_prob):
    if get_fsel_prob is None:
        return None
    return lambda : get_fsel_prob().fsel

@q.timer
def run_psel_from_psel_prob(get_psel_prob):
    if get_psel_prob is None:
        return None
    return lambda : get_psel_prob().psel

# -----------------------------------------------------------------------------

@q.timer_verbose
def run_fsel_prob_sub_sampling(
        job_tag, traj,
        *,
        sub_sampling_rate,
        get_fsel_prob,
        get_f_rand_01,
        get_f_weight,
        ):
    """
    `sub_sampling_rate == 1` implies complete sub-sampling.
    Approximately `sub_sampling_rate` portion of the original selection get selected.
    #
    If `get_f_weight is None` then use `fsel_prob * sub_sampling_rate` as prob to select.
    This is not exactly the same as use `f_weight`!
    """
    assert 1.0 >= sub_sampling_rate >= 0.0
    @q.lazy_call
    @q.timer_verbose
    def get_fsel_prob_sub():
        fname = q.get_fname()
        fsel_prob = get_fsel_prob()
        fsel = fsel_prob.fsel
        f_rand_01 = get_f_rand_01()
        geo = f_rand_01.geo
        sel = fsel[:] >= 0
        f_prob = q.FieldRealD(geo, 1)
        f_prob.set_zero()
        if get_f_weight is not None:
            fsel_rate = get_param(job_tag, "field_selection_fsel_rate")
            f_weight = get_f_weight()
            f_prob @= f_weight
            f_prob *= fsel_rate * sub_sampling_rate
        else:
            f_prob @= fsel_prob
            f_prob *= sub_sampling_rate
        sel_sub = f_prob[:].ravel() >= f_rand_01[:].ravel()
        assert np.all(sel_sub == sel_sub & (fsel[:] >= 0))
        fsel_sub = q.FieldSelection(geo)
        fsel_sub[sel_sub] = fsel[sel_sub]
        fsel_sub.update()
        original_num = q.glb_sum(fsel.n_elems)
        expect_num = original_num * sub_sampling_rate
        actual_num = q.glb_sum(fsel_sub.n_elems)
        q.displayln_info(-1, f"{fname}: sub_sampling_rate = {sub_sampling_rate} ; expect_num = {expect_num} ; actual_num = {actual_num}")
        fsel_prob_sub = q.SelectedFieldRealD(fsel_sub, 1)
        fsel_prob_sub @= f_prob
        return fsel_prob_sub
    return get_fsel_prob_sub

@q.timer_verbose
def run_psel_prob_sub_sampling(
        job_tag, traj,
        *,
        sub_sampling_rate,
        get_psel_prob,
        get_f_rand_01,
        get_f_weight,
        ):
    """
    `sub_sampling_rate == 1` implies complete sub-sampling.
    Approximately `sub_sampling_rate` portion of the original selection get selected.
    #
    If `get_f_weight is None` then use `psel_prob * sub_sampling_rate` as prob to select.
    This is not exactly the same as use `f_weight`!
    This is due to `psel_prob` is always less or equal to 1, `f_weight` does not have any upper limit.
    #
    if get_param(job_tag, "use_simple_psel_sub_sampling", default=False):
        Will simply use the first `sub_sampling_rate * original_num` of the original points.
    """
    assert 1.0 >= sub_sampling_rate >= 0.0
    @q.lazy_call
    @q.timer_verbose
    def get_psel_prob_sub():
        fname = q.get_fname()
        psel_prob = get_psel_prob()
        psel = psel_prob.psel
        original_num = psel.n_points
        total_site = psel.total_site
        use_simple_psel_sub_sampling = get_param(job_tag, "use_simple_psel_sub_sampling", default=False)
        if use_simple_psel_sub_sampling:
            prob_arr = psel_prob[:]
            avg = np.average(prob_arr)
            assert np.linalg.norm(prob_arr - avg) <= 1e-8
            expect_num = round(sub_sampling_rate * original_num)
            xg_arr_sub = psel.xg_arr[:actual_num]
            psel_sub = q.PointsSelection(total_site, xg_arr_sub)
            psel_prob_sub = q.SelectedFieldRealD(psel_sub, 1)
            psel_prob_sub[:] = avg * (expect_num / original_num)
            return psel_prob_sub
        f_rand_01 = get_f_rand_01()
        sp_rand_01 = q.SelectedPointsRealD(psel, 1)
        sp_rand_01 @= f_rand_01
        sp_prob = q.SelectedPointsRealD(psel, 1)
        sp_prob.set_zero()
        if get_f_weight is not None:
            psel_rate = get_param(job_tag, "field_selection_psel_rate")
            f_weight = get_f_weight()
            sp_prob @= f_weight
            sp_prob *= psel_rate * sub_sampling_rate
        else:
            sp_prob @= psel_prob
            sp_prob *= sub_sampling_rate
        sel_sub = sp_prob[:].ravel() >= sp_rand_01[:].ravel()
        xg_arr = psel.xg_arr
        xg_arr_sub = xg_arr[sel_sub]
        psel_sub = q.PointsSelection(total_site, xg_arr_sub)
        expect_num = original_num * sub_sampling_rate
        actual_num = psel_sub.n_points
        q.displayln_info(-1, f"{fname}: sub_sampling_rate = {sub_sampling_rate} ; expect_num = {expect_num} ; actual_num = {actual_num}")
        psel_prob_sub = q.SelectedPointsRealD(psel_sub, 1)
        psel_prob_sub @= sp_prob
        return psel_prob_sub
    return get_psel_prob_sub

# -----------------------------------------------------------------------------

@q.timer(is_timer_fork=True)
def run_psel_split(
        job_tag, traj,
        *,
        get_psel,
        num_piece,
        ):
    """
    Should set `num_piece` be in the form as 2^n.
    return psel_list
    #
    len(psel_list) == num_piece
    """
    assert num_piece >= 1
    fname = q.get_fname()
    path_psel_list = f"{job_tag}/points-selection-split/traj-{traj}/num-piece-{num_piece}"
    @q.lazy_call
    @q.timer(is_timer_fork=True)
    def get_psel_list():
        psel_list = [ q.PointsSelection() for idx in range(num_piece) ]
        load_path = get_load_path(path_psel_list + ".qar")
        assert load_path is not None
        for idx in range(num_piece):
            psel_list[idx].load(get_load_path(f"{path_psel_list}/idx-piece-{idx}.lati"))
        return psel_list
    ret = get_psel_list
    if get_load_path(path_psel_list + "/checkpoint.txt") is not None:
        return ret
    if get_psel is None:
        q.displayln_info(-1, f"{fname}: get_psel is None")
        return None
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    seed = get_job_seed(job_tag)
    rs = q.RngState(f"{seed}-{traj}-{fname}")
    psel = get_psel()
    psel_list = q.psel_split_n_that_increase_separation(psel, num_piece, rs=rs)
    assert len(psel_list) == num_piece
    if 0 == q.get_id_node():
        qar = q.open_qar(get_save_path(path_psel_list + ".qar"), "w")
        for idx in range(num_piece):
            fn = f"idx-piece-{idx}.lati"
            qar.write(fn, "", psel_list[idx].save_str(), skip_if_exist=True)
        qar.write("checkpoint.txt", "", "", skip_if_exist=True)
    q.sync_node()
    q.release_lock()
    return ret

@q.timer(is_timer_fork=True)
def run_fsel_split(
        job_tag, traj,
        *,
        get_fsel,
        num_piece,
        ):
    """
    Should set `num_piece` be in the form as 2^n.
    return psel_list
    #
    len(psel_list) == num_piece
    """
    assert num_piece >= 1
    fname = q.get_fname()
    path_psel_list = f"{job_tag}/field-selection-split/traj-{traj}/num-piece-{num_piece}"
    @q.lazy_call
    @q.timer(is_timer_fork=True)
    def get_psel_list():
        psel_list = [ q.PointsSelection() for idx in range(num_piece) ]
        load_path = get_load_path(path_psel_list + ".qar")
        assert load_path is not None
        for idx in range(num_piece):
            psel_list[idx].load(get_load_path(f"{path_psel_list}/idx-piece-{idx}.lati"))
        return psel_list
    ret = get_psel_list
    if get_load_path(path_psel_list + "/checkpoint.txt") is not None:
        return ret
    if get_fsel is None:
        q.displayln_info(-1, f"{fname}: get_fsel is None")
        return None
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}"):
        return None
    seed = get_job_seed(job_tag)
    rs = q.RngState(f"{seed}-{traj}-{fname}")
    psel = get_fsel().to_psel()
    psel_list = q.psel_split_n_that_increase_separation(psel, num_piece, rs=rs)
    assert len(psel_list) == num_piece
    if 0 == q.get_id_node():
        qar = q.open_qar(get_save_path(path_psel_list + ".qar"), "w")
        for idx in range(num_piece):
            fn = f"idx-piece-{idx}.lati"
            qar.write(fn, "", psel_list[idx].save_str(), skip_if_exist=True)
        qar.write("checkpoint.txt", "", "", skip_if_exist=True)
    q.sync_node()
    q.release_lock()
    return ret

# -----------------------------------------------------------------------------

@q.timer_verbose
def save_prop_wsrc_sparse(job_tag, traj, *, load_prop, tslice, inv_type, inv_acc, sfw, qar_sp, psel, fsel):
    """
    `load_prop()` is the wall source propagator, as if generated in a Coulomb gauge fixed gauge configuration.
    `tslice` is the integer indicate the source time slice.
    `inv_type` is `0` (light quark) or `1` (strange quark) or `2` (charm quark)
    `inv_acc` for wall source typically be integer `1` (sloppy accuracy) or `2` (exact accuracy).
    """
    fname = q.get_fname()
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in sfw:
        if 0 == q.get_id_node():
            assert f"{tag}.lat" in qar_sp
            assert f"{tag} ; wsnk.lat" in qar_sp
        q.displayln_info(0, f"{fname}: tag='{tag}' of '{job_tag}/{traj}' already saved. Skipping.")
        return
    prop = load_prop()
    if prop is None:
        return
    q.displayln_info(-1, f"{fname}: write wsrc prop tag='{tag}' of '{job_tag}/{traj}'")
    s_prop = q.SelProp(fsel)
    ps_prop = q.PselProp(psel)
    q.set_zero(s_prop)
    q.set_zero(ps_prop)
    s_prop @= prop
    ps_prop @= prop
    ps_prop_ws = prop.glb_sum_tslice()
    qar_sp.write(f"{tag}.lat", "", ps_prop.save_str(), skip_if_exist=True)
    qar_sp.write(f"{tag} ; wsnk.lat", "", ps_prop_ws.save_str(), skip_if_exist=True)
    s_prop.save_float_from_double(sfw, tag, skip_if_exist=True)
    qar_sp.flush()
    sfw.flush()

@q.timer(is_timer_fork=True)
def run_prop_wsrc_sparse(
        job_tag, traj,
        *,
        inv_type, get_gf, get_gt, get_eig, get_psel, get_fsel, get_wi,
        ):
    fname = q.get_fname()
    if None in [ get_gt, get_psel, get_fsel, get_wi, ]:
        return
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    path_f = f"{job_tag}/prop-wsrc-full-{inv_type_name}/traj-{traj}/geon-info.txt"
    path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-wsrc-{inv_type_name}/traj-{traj}"
    is_performing_inversion = get_param(job_tag, "is_performing_inversion_if_no_full_prop_available", default=False)
    if get_load_path(path_f) is None:
        if not is_performing_inversion:
            q.displayln_info(f"WARNING: {fname}: {job_tag} {traj} {inv_type_name} full prop is not available yet.")
            return
    if get_load_path(path_s + "/geon-info.txt") is not None:
        assert get_load_path(path_sp + "/checkpoint.txt") is not None
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    gt = get_gt()
    fsel = get_fsel()
    psel = get_psel()
    wi = get_wi()
    if get_load_path(path_f) is not None:
        sfr = q.open_fields(get_load_path(path_f), "r")
        available_tags = sfr.list()
        q.displayln_info(0, f"available_tags={available_tags}")
    else:
        assert is_performing_inversion
        sfr = None
        if None in [ get_gf, get_eig, ]:
            assert False
        gf = get_gf()
        eig = get_eig()
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    qar_sp = q.open_qar_info(get_save_path(path_sp + ".qar"), "a")
    for idx, tslice, inv_type_wi, inv_acc in wi:
        if inv_type_wi != inv_type:
            continue
        def load_prop():
            tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
            prop = q.Prop(geo)
            if sfr is not None:
                if tag not in sfr:
                    raise Exception(f"{tag} not in {sfr.list()}")
                prop.load_double(sfr, tag)
            else:
                prop = compute_prop_wsrc_1(
                        job_tag, traj,
                        gf=gf, gt=gt, eig=eig,
                        idx=idx, tslice=tslice,
                        inv_type=inv_type, inv_acc=inv_acc,
                        )
            return prop
        save_prop_wsrc_sparse(
                job_tag, traj,
                load_prop=load_prop,
                tslice=tslice, inv_type=inv_type, inv_acc=inv_acc,
                sfw=sfw, qar_sp=qar_sp,
                psel=psel, fsel=fsel,
                )
    sfw.close()
    qar_sp.write("checkpoint.txt", "", "", skip_if_exist=True)
    qar_sp.flush()
    qar_sp.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.release_lock()
    return [ f"{fname} {job_tag} {traj} {inv_type} done", ]

# -----------------------------------------------------------------------------

@q.timer
def calc_hvp_sum_tslice(chvp_16):
    """
    return ld_hvp_ts
    #
    ld_hvp_ts[t_dir, t, mu, nu]
    #
    t_dir in [ "x", "y", "z", "t", ]
    t in range(t_size)
    mu in range(4)
    nu in range(4)
    #
    `mu` for sink and `nu` for source
    t_size = max(total_site)
    #
    arr_for_t_dir.shape == (t_size, 4, 4,)
    #
    chvp_16 is a hvp field (from q.contract_chvp_16)
    """
    total_site = chvp_16.total_site
    t_size = max(total_site)
    ld_hvp_ts = q.mk_lat_data([
        [ "t_dir", 4, [ "x", "y", "z", "t", ], ],
        [ "t", t_size, ],
        [ "mu", 4, ],
        [ "nu", 4, ],
        ])
    ld_hvp_ts.set_zero()
    ld_arr = ld_hvp_ts[:]
    assert ld_arr.shape == (4, t_size, 4, 4,)
    assert ld_arr.dtype == np.complex128
    for t_dir in range(4):
        chvp_16_ts = chvp_16.glb_sum_tslice(t_dir=t_dir)
        arr = chvp_16_ts.to_numpy()
        t_size = arr.shape[0]
        ld_arr[t_dir, :t_size] = arr.reshape(t_size, 4, 4)
    return ld_hvp_ts

@q.timer
def compute_prop_psrc_hvp_contract(
        job_tag, traj, xg_src, inv_type, inv_acc,
        *,
        prop, tag, sfw_hvp, qar_hvp_ts):
    """
    # chvp_16.get_elem(x, mu * 4 + nu) is complex
    # (1) mu is the sink polarization and nu is the src polarization
    # (2) hvp field is simply the trace of the products of gamma matrix and propagators.
    #     It does not include the any minus sign (e.g. The minus sign due to the loop).
    """
    if (sfw_hvp is None) and (qar_hvp_ts is None):
        return
    assert isinstance(xg_src, q.Coordinate)
    chvp_16 = q.contract_chvp_16(prop, prop)
    ld_hvp_ts = calc_hvp_sum_tslice(chvp_16)
    if qar_hvp_ts is not None:
        qar_hvp_ts.write(f"{tag}.lat", "", ld_hvp_ts.save_str(), skip_if_exist=True)
        qar_hvp_ts.flush()
    if sfw_hvp is not None:
        chvp_16.save_float_from_double(sfw_hvp, tag, skip_if_exist=True)
        sfw_hvp.flush()

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_prop_2(inv, src, *, tag, sfw, qar_sp, psel, fsel,
                   f_rand_01, fsel_psrc_prop_norm_threshold, gt):
    sol = inv * src
    sp_sol = q.PselProp(psel)
    sp_sol @= sol
    qar_sp.write(f"{tag}.lat", "", sp_sol.save_str(), skip_if_exist=True)
    sol_gt = gt * sol
    sol_ws = sol_gt.glb_sum_tslice()
    qar_sp.write(f"{tag} ; wsnk.lat", "", sol_ws.save_str(), skip_if_exist=True)
    sol_ps_sel_prob = q.qnorm_field(sol)
    sol_ps_sel_prob *= 1.0 / fsel_psrc_prop_norm_threshold
    sol_ps_sel_prob[:] = np.minimum(1.0, sol_ps_sel_prob[:])
    ps_sel = f_rand_01[:, 0] <= sol_ps_sel_prob[:, 0]
    fsel_ps = q.FieldSelection(fsel.geo)
    fsel_ps[ps_sel] = 0
    fsel_ps.update()
    fsel_combine = fsel_ps.copy()
    fsel_combine.add_fsel(fsel)
    num_fsel = q.glb_sum(fsel.n_elems)
    num_fsel_ps = q.glb_sum(fsel_ps.n_elems)
    num_fsel_combine = q.glb_sum(fsel_combine.n_elems)
    q.displayln_info(0, f"compute_prop_psrc: tag='{tag}' ; num_fsel={num_fsel} ; num_fsel_ps={num_fsel_ps} ; num_fsel_combine={num_fsel_combine}")
    s_sol_ps_sel_prob = q.SelectedFieldRealD(fsel_ps)
    s_sol_ps_sel_prob @= sol_ps_sel_prob
    s_sol_ps_sel_prob.save_double(sfw, f"{tag} ; fsel-prob-psrc-prop", skip_if_exist=True)
    s_sol = q.SelProp(fsel_combine)
    s_sol @= sol
    s_sol.save_float_from_double(sfw, tag, skip_if_exist=True)
    qar_sp.flush()
    sfw.flush()
    return sol

def mk_psrc_tag(xg, inv_type, inv_acc):
    xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
    tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
    return tag

@q.timer
def compute_prop_psrc(job_tag, traj, xg_src, inv_type, inv_acc, *,
        idx, gf, gt, sfw, qar_sp, psel, fsel, f_rand_01, sfw_hvp, qar_hvp_ts,
        eig):
    assert isinstance(xg_src, q.Coordinate)
    tag = mk_psrc_tag(xg_src, inv_type, inv_acc)
    if (tag in sfw) and (sfw_hvp is None or (tag in sfw_hvp)):
        assert f"{tag} ; fsel-prob-psrc-prop" in sfw
        if 0 == q.get_id_node():
            assert f"{tag}.lat" in qar_sp
            assert f"{tag} ; wsnk.lat" in qar_sp
        if qar_hvp_ts is not None:
            assert qar_hvp_ts.has_regular_file(f"{tag}.lat")
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_psrc: {job_tag} idx={idx} tag='{tag}'")
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig=eig)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    fsel_psrc_prop_norm_threshold = get_param(job_tag, "field_selection_fsel_psrc_prop_norm_threshold")
    geo = q.Geometry(total_site)
    src = q.mk_point_src(geo, xg_src)
    prop = compute_prop_2(
            inv, src, tag=tag, sfw=sfw, qar_sp=qar_sp,
            psel=psel, fsel=fsel,
            f_rand_01=f_rand_01,
            fsel_psrc_prop_norm_threshold=fsel_psrc_prop_norm_threshold,
            gt=gt)
    compute_prop_psrc_hvp_contract(
            job_tag, traj, xg_src, inv_type, inv_acc,
            prop=prop, tag=tag, sfw_hvp=sfw_hvp, qar_hvp_ts=qar_hvp_ts)

@q.timer_verbose
def compute_prop_psrc_all(job_tag, traj, *,
                          inv_type, gf, gt, psel, fsel, f_rand_01, eig):
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    path_s = f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}"
    path_s_hvp = f"{job_tag}/hvp-psrc-{inv_type_name}/traj-{traj}"
    path_hvp_ts = f"{job_tag}/hvp-sum-tslice-psrc-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-psrc-{inv_type_name}/traj-{traj}"
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    is_saving_hvp = get_param(job_tag, "run_prop_psrc", "is_saving_hvp", default=True)
    is_saving_hvp_ts = get_param(job_tag, "run_prop_psrc", "is_saving_hvp_ts", default=True)
    if is_saving_hvp:
        sfw_hvp = q.open_fields(get_save_path(path_s_hvp + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    else:
        sfw_hvp = None
    qar_sp = q.open_qar_info(get_save_path(path_sp + ".qar"), "a")
    if is_saving_hvp_ts:
        qar_hvp_ts = q.open_qar_info(get_save_path(path_hvp_ts + ".qar"), "a")
    else:
        qar_hvp_ts = None
    def comp(idx, xg_src, inv_acc):
        compute_prop_psrc(job_tag, traj, xg_src, inv_type, inv_acc,
                idx=idx, gf=gf, gt=gt, sfw=sfw, qar_sp=qar_sp,
                psel=psel, fsel=fsel, f_rand_01=f_rand_01,
                sfw_hvp=sfw_hvp, qar_hvp_ts=qar_hvp_ts,
                eig=eig)
    prob1 = get_param(job_tag, "prob_acc_1_psrc")
    prob2 = get_param(job_tag, "prob_acc_2_psrc")
    seed = get_job_seed(job_tag)
    rs = q.RngState(f"seed {seed} {traj}").split(f"compute_prop_psrc_all(ama)")
    for idx, xg_src in enumerate(psel):
        r = rs.split(f"{xg_src.to_tuple()}").u_rand_gen()
        assert 0 <= r and r <= 1
        comp(idx, xg_src, inv_acc=0)
        if r <= prob1:
            comp(idx, xg_src, inv_acc=1)
        if r <= prob2:
            comp(idx, xg_src, inv_acc=2)
    sfw.close()
    if sfw_hvp is not None:
        sfw_hvp.close()
    qar_sp.write("checkpoint.txt", "", "", skip_if_exist=True)
    qar_sp.flush()
    qar_sp.close()
    if qar_hvp_ts is not None:
        qar_hvp_ts.write("checkpoint.txt", "", "", skip_if_exist=True)
        qar_hvp_ts.flush()
        qar_hvp_ts.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    if sfw_hvp is not None:
        q.qrename_info(get_save_path(path_s_hvp + ".acc"), get_save_path(path_s_hvp))

@q.timer(is_timer_fork=True)
def run_prop_psrc(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt, get_psel, get_fsel, get_f_rand_01):
    fname = q.get_fname()
    if None in [ get_gf, get_gt, get_psel, get_fsel, get_f_rand_01, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    if get_load_path(f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        psel = get_psel()
        fsel = get_fsel()
        f_rand_01 = get_f_rand_01()
        assert fsel.is_containing(psel)
        compute_prop_psrc_all(job_tag, traj,
                              inv_type=inv_type, gf=gf, gt=gt,
                              psel=psel, fsel=fsel,
                              f_rand_01=f_rand_01,
                              eig=eig)
        q.release_lock()
        return [ f"{fname} {job_tag} {traj} {inv_type} done", ]

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_hvp_average(job_tag, traj, *, inv_type, psel_prob, data_path, geo):
    fname = q.get_fname()
    psel = psel_prob.psel
    psel_prob_arr = psel_prob[:].ravel()
    hvp_average = q.FieldComplexD(geo, 16)
    hvp_average.set_zero()
    rel_acc_list = [ 0, 1, 2, ]
    prob_list = [ get_param(job_tag, f"prob_acc_{inv_acc}_psrc") for inv_acc in rel_acc_list ]
    sfr = q.open_fields(data_path, "r")
    tags = sfr.list()
    for p_idx in range(len(psel)):
        xg_src = q.Coordinate(psel[p_idx])
        prob_src = psel_prob_arr[p_idx]
        val_list = []
        for inv_acc in rel_acc_list:
            tag = mk_psrc_tag(xg_src, inv_type, inv_acc)
            if tag not in tags:
                val_list.append(None)
            else:
                chvp_16 = q.FieldComplexD(geo, 16)
                chvp_16.load_double_from_float(sfr, tag)
                val_list.append(chvp_16)
        assert val_list[0] is not None
        ama_val = q.mk_ama_val(val_list[0], xg_src.to_tuple(), val_list, rel_acc_list, prob_list)
        hvp = q.ama_extract(ama_val).shift(-xg_src)
        hvp *= 1 / prob_src
        hvp_average += hvp
    sfr.close()
    hvp_average *= 1 / geo.total_volume
    return hvp_average

@q.timer(is_timer_fork=True)
def run_hvp_average(job_tag, traj, *, inv_type, get_psel_prob):
    """
    return get_hvp_average()
    save hvp_average.field in single precision.
    hvp_average.get_elem(x, mu * 4 + nu) is complex
    #
    (1) mu is the sink polarization and nu is the src polarization
    (2) hvp field is simply the trace of the products of gamma matrix and propagators.
        It does not include the any minus sign (e.g. The minus sign due to the loop).
    """
    fname = q.get_fname()
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    fn = f"{job_tag}/hvp-average/traj-{traj}/hvp_average_{inv_type_name}.field"
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    @q.lazy_call
    def load():
        hvp_average = q.FieldComplexD(geo, 16)
        hvp_average.load_double_from_float(get_load_path(fn))
        return hvp_average
    if get_load_path(fn) is not None:
        return load
    data_fn = f"{job_tag}/hvp-psrc-{inv_type_name}/traj-{traj}/geon-info.txt"
    data_path = get_load_path(data_fn)
    if data_path is None:
        q.displayln_info(f"{fname}: '{data_fn}' not available.")
        return
    q.check_stop()
    q.check_time_limit()
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        return
    psel_prob = get_psel_prob()
    hvp_average = compute_hvp_average(job_tag, traj, inv_type=inv_type, psel_prob=psel_prob, geo=geo, data_path=data_path)
    hvp_average.save_float_from_double(get_save_path(fn))
    q.displayln_info(f"{fname}: {job_tag} {traj} {inv_type} num_fields={len(psel_prob.psel)}")
    q.release_lock()
    return load

# -----------------------------------------------------------------------------

@q.timer(is_timer_fork=True)
def run_field_rand_u1_dict(
        job_tag, traj,
        ):
    """
    return get_field_rand_u1_dict
    #
    get_field_rand_u1_dict()["fsel-src"] => q.FieldComplexD
    get_field_rand_u1_dict()["fsel-src-dag"] => q.FieldComplexD
    get_field_rand_u1_dict()["psel-src"] => q.FieldComplexD
    get_field_rand_u1_dict()["psel-src-dag"] => q.FieldComplexD
    """
    fname = q.get_fname()
    path = f"{job_tag}/field-rand-u1/traj-{traj}"
    @q.lazy_call
    @q.timer_verbose
    def get_field_rand_u1_dict():
        d = dict()
        for psel_list_type in [ "fsel", "psel" ]:
            fn = get_load_path(f"{path}/{psel_list_type}-src.field")
            fn_dag = get_load_path(f"{path}/{psel_list_type}-src-dag.field")
            assert fn is not None
            assert fn_dag is not None
            f = q.FieldComplexD()
            f_dag = q.FieldComplexD()
            f.load_double(fn)
            f_dag.load_double(fn_dag)
            d[f"{psel_list_type}-src"] = f
            d[f"{psel_list_type}-src-dag"] = f_dag
        return d
    ret = get_field_rand_u1_dict
    if get_load_path(f"{path}/checkpoint.txt"):
        return ret
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    seed = get_job_seed(job_tag)
    rs_rand_u1 = q.RngState(f"seed {seed} {traj}").split(f"compute_prop_rand_sparse_u1_src(rand_u1)")
    for psel_list_type in [ "fsel", "psel" ]:
        fu1 = q.mk_rand_vol_u1(geo, rs_rand_u1.split(f"{psel_list_type}"))
        fu1_dag = q.FieldComplexD(geo, 1)
        fu1_dag[:] = fu1[:].conj()
        fu1.save_double(get_save_path(f"{path}/{psel_list_type}-src.field"))
        fu1_dag.save_double(get_save_path(f"{path}/{psel_list_type}-src-dag.field"))
        if is_test():
            q.json_results_append(f"{fname}: {psel_list_type} fu1", q.get_data_sig_arr(fu1, q.RngState(), 4), 1e-15)
            q.json_results_append(f"{fname}: {psel_list_type} fu1_dag", q.get_data_sig_arr(fu1_dag, q.RngState(), 4), 1e-15)
    q.qtouch_info(get_save_path(f"{path}/checkpoint.txt"), "")
    return ret

@q.timer(is_timer_fork=True)
def run_prop_sparse_rand_u1_src(
        job_tag, traj,
        *,
        inv_type,
        get_gf,
        get_psel,
        get_fsel,
        get_field_rand_u1_dict,
        get_psel_list=None,
        get_fsel_psel_list=None,
        get_eig=None,
        ):
    """
    fsel should contain psel
    Should set one (and only one) of the `get_psel_list` and `get_fsel_psel_list`.
    """
    fname = q.get_fname()
    if None in [ get_gf, get_psel, get_fsel, get_field_rand_u1_dict, ]:
        return
    if (get_psel_list is None) and (get_fsel_psel_list is None):
        return None
    if get_eig is None:
        if inv_type == 0:
            return None
        else:
            get_eig = lambda: None
    if get_psel_list is not None:
        assert get_fsel_psel_list is None
        psel_list_type = "psel"
    elif get_fsel_psel_list is not None:
        psel_list_type = "fsel"
    else:
        assert False
    quark_flavor_list = get_param(job_tag, "quark_flavor_list")
    quark_flavor = quark_flavor_list[inv_type]
    path_s = f"{job_tag}/prop-rand-u1-{psel_list_type}-sparse-{quark_flavor}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-rand-u1-{psel_list_type}-sparse-{quark_flavor}/traj-{traj}"
    if get_load_path(path_s + "/geon-info.txt") is not None:
        assert get_load_path(path_sp + "/checkpoint.txt") is not None
        return
    if not q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{psel_list_type}-{quark_flavor}"):
        return
    gf = get_gf()
    geo = gf.geo
    fsel = get_fsel()
    psel = get_psel()
    field_rand_u1_dict = get_field_rand_u1_dict()
    eig = get_eig()
    if psel_list_type == "psel":
        assert get_psel_list is not None
        assert get_fsel_psel_list is None
        psel_list = get_psel_list()
        fu1 = field_rand_u1_dict["psel-src"]
        fu1_dag = field_rand_u1_dict["psel-src-dag"]
    elif psel_list_type == "fsel":
        assert get_psel_list is None
        assert get_fsel_psel_list is not None
        psel_list = get_fsel_psel_list()
        fu1 = field_rand_u1_dict["fsel-src"]
        fu1_dag = field_rand_u1_dict["fsel-src-dag"]
    else:
        assert False
    total_site = geo.total_site
    gt = None
    prob_acc_1_rand_sparse_u1 = get_param(job_tag, "prob_acc_1_rand_u1_sparse")
    prob_acc_2_rand_sparse_u1 = get_param(job_tag, "prob_acc_2_rand_u1_sparse")
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    qar_sp = q.open_qar_info(get_save_path(path_sp + ".qar"), "a")
    fu1_from_sp = q.FieldComplexD(geo, 1)
    seed = get_job_seed(job_tag)
    rs_ama = q.RngState(f"seed {seed} {traj}").split(f"compute_prop_rand_u1(ama)")
    @q.timer
    def compute_and_save(idx_psel, is_dagger, inv_acc):
        nonlocal fu1_from_sp
        assert 0 <= idx_psel < len(psel_list)
        assert is_dagger in [ 0, 1, ]
        tag = f"idx_psel={idx_psel} ; is_dagger={is_dagger} ; type={inv_type} ; accuracy={inv_acc}"
        if tag in sfw:
            if q.get_id_node() == 0:
                assert f"{tag} ; psel_rand_src.lati" in qar_sp
                assert f"{tag} ; psel_rand_src ; fu1.lat" in qar_sp
                assert f"{tag} ; psel_rand_src ; prop.lat" in qar_sp
                assert f"{tag}.lat" in qar_sp
            return
        if q.get_id_node() == 0:
            if f"{tag}.lat" in qar_sp:
                q.displayln_info(-1, f"WARNING: {fname}: '{tag}.lat' already exist in {qar_sp.path()}")
        if f"{tag} ; fu1" in sfw:
            q.displayln_info(-1, f"WARNING: {fname}: '{tag} ; fu1' already exist in {sfw.path()}")
        q.check_stop()
        q.check_time_limit()
        psel_rand_src = psel_list[idx_psel]
        sp_fu1 = q.SelectedPointsComplexD(psel_rand_src, 1)
        sp_fu1.set_zero()
        if is_dagger == 0:
            sp_fu1 @= fu1
        elif is_dagger == 1:
            sp_fu1 @= fu1_dag
        else:
            assert False
        fu1_from_sp.set_zero()
        fu1_from_sp @= sp_fu1
        prop_src = q.mk_rand_vol_u1_src(fu1_from_sp)
        inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)
        prop_sol = inv * prop_src
        s_prop = q.SelProp(fsel)
        ps_prop = q.PselProp(psel)
        ps_prop_prs = q.PselProp(psel_rand_src)
        s_prop @= prop_sol
        ps_prop @= prop_sol
        ps_prop_prs @= prop_sol
        qar_sp.write(f"{tag} ; psel_rand_src.lati", "", psel_rand_src.save_str(), skip_if_exist=True)
        qar_sp.write(f"{tag} ; psel_rand_src ; fu1.lat", "", sp_fu1.save_str(), skip_if_exist=True)
        qar_sp.write(f"{tag} ; psel_rand_src ; prop.lat", "", ps_prop_prs.save_str(), skip_if_exist=True)
        qar_sp.write(f"{tag}.lat", "", ps_prop.save_str(), skip_if_exist=True)
        qar_sp.flush()
        s_prop.save_float_from_double(sfw, tag, skip_if_exist=True)
        sfw.flush()
    for idx_psel in range(len(psel_list)):
        r = rs_ama.split(str(idx_psel)).u_rand_gen()
        assert 0 <= r and r <= 1
        for is_dagger in [ 0, 1, ]:
            inv_acc = 0
            compute_and_save(idx_psel, is_dagger, inv_acc)
            if r <= prob_acc_1_rand_sparse_u1:
                inv_acc = 1
                compute_and_save(idx_psel, is_dagger, inv_acc)
            if r <= prob_acc_2_rand_sparse_u1:
                inv_acc = 2
                compute_and_save(idx_psel, is_dagger, inv_acc)
    sfw.close()
    qar_sp.write("checkpoint.txt", "", "", skip_if_exist=True)
    qar_sp.flush()
    qar_sp.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.clean_cache(q.cache_inv)
    if is_test():
        q.json_results_append(f"{fname} {job_tag} {traj} inv_type={inv_type} psel_list_type={psel_list_type}")
        sfr = q.open_fields(get_load_path(path_s), "r")
        sfr_list = sfr.list()
        q.json_results_append(f"{fname} sfr.list()={sfr_list}")
        for tag in sfr_list:
            s_prop = q.SelProp()
            s_prop.load_double_from_float(sfr, tag)
            q.json_results_append(f"{fname} {tag}", q.get_data_sig_arr(s_prop, q.RngState(), 2), 1e-4)
        sfr.close()
        qar_sp = q.open_qar_info(get_load_path(path_sp + ".qar"), "r")
        qar_sp_list = qar_sp.list()
        qar_sp_list = q.get_comm().bcast(qar_sp_list)
        q.json_results_append(f"{fname} qar_sp.list()={qar_sp_list}")
        for tag in qar_sp_list:
            if tag.endswith(".lati"):
                ld = q.load_lat_data_int(get_load_path(f"{path_sp}/{tag}"))
            elif tag.endswith(".lat"):
                ld = q.load_lat_data(get_load_path(f"{path_sp}/{tag}"))
            elif tag == "checkpoint.txt":
                continue
            else:
                assert False
            q.json_results_append(f"{fname} {tag}", q.get_data_sig_arr(ld, q.RngState(), 2), 1e-4)
        qar_sp.close()
    q.release_lock()
    return [ f"{fname} {job_tag} {traj} {inv_type} done", ]

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_prop_rand_u1_type_acc(*, sfw, job_tag, traj, gf, eig, fsel, idx_rand_u1, inv_type, inv_acc):
    """
    same rand source for different inv_type
    """
    tag = f"idx_rand_u1={idx_rand_u1} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in sfw:
        return
    q.check_stop()
    q.check_time_limit()
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig=eig)
    seed = get_job_seed(job_tag)
    rs = q.RngState(f"seed {seed} {traj}").split(f"compute_prop_rand_u1(rand_u1)").split(str(idx_rand_u1))
    s_prop = q.mk_rand_u1_prop(inv, fsel, rs)
    s_prop.save_float_from_double(sfw, tag, skip_if_exist=True)
    sfw.flush()
    return s_prop

@q.timer_verbose
def compute_prop_rand_u1(*, job_tag, traj, inv_type, gf, path_s, fsel, eig=None):
    """
    use fsel instead of fselc
    """
    n_rand_u1_fsel = get_param(job_tag, "n_rand_u1_fsel")
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    def comp(idx_rand_u1, inv_acc):
        compute_prop_rand_u1_type_acc(
                sfw=sfw,
                job_tag=job_tag, traj=traj,
                gf=gf, eig=eig, fsel=fsel,
                idx_rand_u1=idx_rand_u1,
                inv_type=inv_type, inv_acc=inv_acc)
    seed = get_job_seed(job_tag)
    rs = q.RngState(f"seed {seed} {traj}").split(f"compute_prop_rand_u1(ama)")
    prob1 = get_param(job_tag, "prob_acc_1_rand_u1")
    prob2 = get_param(job_tag, "prob_acc_2_rand_u1")
    for idx_rand_u1 in range(n_rand_u1_fsel):
        r = rs.split(str(idx_rand_u1)).u_rand_gen()
        inv_acc = 0
        assert 0 <= r and r <= 1
        comp(idx_rand_u1, inv_acc)
        inv_acc = 1
        if r <= prob1:
            comp(idx_rand_u1, inv_acc)
        inv_acc = 2
        if r <= prob2:
            comp(idx_rand_u1, inv_acc)
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

@q.timer(is_timer_fork=True)
def run_prop_rand_u1(job_tag, traj, *, inv_type, get_gf, get_fsel, get_eig=None):
    fname = q.get_fname()
    if None in [ get_gf, get_fsel, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_name_list = get_param(job_tag, "quark_flavor_list", default=[ "light", "strange", "charm", ])
    inv_type_name = inv_type_name_list[inv_type]
    path_s = f"{job_tag}/prop-rand-u1-{inv_type_name}/traj-{traj}"
    if get_load_path(path_s + "/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        gf = get_gf()
        fsel = get_fsel()
        eig = get_eig()
        compute_prop_rand_u1(
                job_tag=job_tag, traj=traj,
                inv_type=inv_type,
                gf=gf,
                path_s=path_s,
                fsel=fsel,
                eig=eig)
        q.release_lock()
        return [ f"{fname} {job_tag} {traj} {inv_type} done", ]

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_prop_3(inv, src_smear, *, tag, sfw, qar_sp, sfw_m, psel, fsel, gt, psel_smear, fsel_smear_median, smear):
    sol = inv * src_smear
    sp_sol = q.PselProp(psel)
    sp_sol @= sol
    qar_sp.write(f"{tag}.lat", "", sp_sol.save_str(), skip_if_exist=True)
    sol_gt = gt * sol
    sol_ws = sol_gt.glb_sum_tslice()
    qar_sp.write(f"{tag} ; wsnk.lat", "", sol_ws.save_str(), skip_if_exist=True)
    sol_smear = smear(sol)
    sol_smear_psel = q.PselProp(psel_smear)
    sol_smear_psel @= sol_smear
    qar_sp.write(f"{tag} ; smear-snk.lat", "", sol_smear_psel.save_str(), skip_if_exist=True)
    sm_sol = q.SelProp(fsel_smear_median)
    sm_sol @= sol_smear
    sm_sol.save_float_from_double(sfw_m, f"{tag} ; smear-snk", skip_if_exist=True)
    s_sol = q.SelProp(fsel)
    s_sol @= sol
    qar_sp.flush()
    sfw_m.flush()
    s_sol.save_float_from_double(sfw, tag, skip_if_exist=True)
    sfw.flush()
    return sol

@q.timer
def compute_prop_smear(job_tag, xg_src, inv_type, inv_acc, *,
        idx, gf, gt, sfw, qar_sp, sfw_m, psel, fsel, psel_smear, fsel_smear_median, gf_ape, eig):
    xg = xg_src
    xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
    tag = f"smear ; xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in sfw:
        if 0 == q.get_id_node():
            assert f"{tag}.lat" in qar_sp
            assert f"{tag} ; wsnk.lat" in qar_sp
            assert f"{tag} ; smear-snk.lat" in qar_sp
        assert f"{tag} ; smear-snk" in sfw_m
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_smear: {job_tag} idx={idx} tag='{tag}'")
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig=eig)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    coef = get_param(job_tag, "prop_smear_coef")
    step = get_param(job_tag, "prop_smear_step")
    def smear(src):
        return q.prop_smear(src, gf_ape, coef, step)
    src = smear(q.mk_point_src(geo, xg_src))
    prop = compute_prop_3(
            inv, src, tag=tag,
            sfw=sfw, qar_sp=qar_sp, sfw_m=sfw_m,
            psel=psel, fsel=fsel, gt=gt,
            psel_smear=psel_smear, fsel_smear_median=fsel_smear_median,
            smear=smear,
            )

@q.timer_verbose
def compute_prop_smear_all(
        job_tag, traj, *,
        inv_type, gf, gt, gf_ape, eig,
        psel, fsel, psel_smear, fsel_smear_median,
        ):
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    path_s = f"{job_tag}/prop-smear-{inv_type_name}/traj-{traj}"
    path_sm = f"{job_tag}/psel_smear_median-prop-smear-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-smear-{inv_type_name}/traj-{traj}"
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    sfw_m = q.open_fields(get_save_path(path_sm + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
    qar_sp = q.open_qar_info(get_save_path(path_sp + ".qar"), "a")
    def comp(idx, xg_src, inv_acc):
        compute_prop_smear(
                job_tag, xg_src, inv_type, inv_acc,
                idx=idx,
                gf=gf, gt=gt, gf_ape=gf_ape, eig=eig,
                sfw=sfw, qar_sp=qar_sp, sfw_m=sfw_m,
                psel=psel, fsel=fsel,
                psel_smear=psel_smear,
                fsel_smear_median=fsel_smear_median,
                )
    prob1 = get_param(job_tag, "prob_acc_1_smear")
    prob2 = get_param(job_tag, "prob_acc_2_smear")
    seed = get_job_seed(job_tag)
    rs = q.RngState(f"seed {seed} {traj}").split(f"compute_prop_smear_all(ama)")
    for idx, xg_src in enumerate(psel_smear):
        r = rs.split(f"{xg_src.to_tuple()}").u_rand_gen()
        assert 0 <= r and r <= 1
        comp(idx, xg_src, inv_acc=0)
        if r <= prob1:
            comp(idx, xg_src, inv_acc=1)
        if r <= prob2:
            comp(idx, xg_src, inv_acc=2)
    sfw.close()
    qar_sp.write("checkpoint.txt", "", "", skip_if_exist=True)
    qar_sp.flush()
    qar_sp.close()
    q.qrename_info(get_save_path(path_sm + ".acc"), get_save_path(path_sm))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

@q.timer(is_timer_fork=True)
def run_prop_smear(job_tag, traj, *, inv_type, get_gf, get_gf_ape, get_eig, get_gt, get_psel, get_fsel, get_psel_smear, get_psel_smear_median):
    fname = q.get_fname()
    if None in [ get_gf, get_gt, get_gf_ape, get_psel, get_fsel, get_psel_smear, get_psel_smear_median, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_name_list = [ "light", "strange", ]
    inv_type_name = inv_type_name_list[inv_type]
    if get_load_path(f"{job_tag}/prop-smear-{inv_type_name}/traj-{traj}/geon-info.txt") is not None:
        assert get_load_path(f"{job_tag}/psel-prop-smear-{inv_type_name}/traj-{traj}/checkpoint.txt") is not None
        assert get_load_path(f"{job_tag}/psel_smear_median-prop-smear-{inv_type_name}/traj-{traj}/geon-info.txt") is not None
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-{fname}-{inv_type_name}"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        psel = get_psel()
        fsel = get_fsel()
        assert fsel.is_containing(psel)
        psel_smear = get_psel_smear()
        psel_smear_median = get_psel_smear_median()
        fsel_smear_median = q.FieldSelection(psel_smear_median)
        assert psel_smear_median.is_containing(psel_smear)
        gf_ape = get_gf_ape()
        compute_prop_smear_all(
                job_tag, traj,
                inv_type=inv_type, gf=gf, gt=gt, gf_ape=gf_ape, eig=eig,
                psel=psel, fsel=fsel, psel_smear=psel_smear, fsel_smear_median=fsel_smear_median,
                )
        q.release_lock()
        return [ f"{fname} {job_tag} {traj} {inv_type} done", ]

# -----------------------------------------------------------------------------

# @q.timer_verbose
# def compute_prop_1(inv, src, *, tag, sfw, path_sp, psel, fsel):
#     fn_sp = os.path.join(path_sp, f"{tag}.lat")
#     fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
#     sol = inv * src
#     sp_sol = q.PselProp(psel)
#     sp_sol @= sol
#     sp_sol.save(get_save_path(fn_sp))
#     sol_ws = sol.glb_sum_tslice()
#     sol_ws.save(get_save_path(fn_spw))
#     s_sol = q.SelProp(fsel)
#     s_sol @= sol
#     s_sol.save_float_from_double(sfw, tag)
#     sfw.flush()
#     return sol
#
# @q.timer
# def compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc, *,
#         idx, sfw, path_sp, psel, fsel, eig, finished_tags):
#     tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
#     if tag in finished_tags:
#         return None
#     q.check_stop()
#     q.check_time_limit()
#     q.displayln_info(f"compute_prop_wsrc: idx={idx} tslice={tslice}", job_tag, inv_type, inv_acc)
#     inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt=gt, eig=eig)
#     total_site = q.Coordinate(get_param(job_tag, "total_site"))
#     geo = q.Geometry(total_site)
#     src = q.mk_wall_src(geo, tslice)
#     prop = compute_prop_1(inv, src, tag=tag, sfw=sfw, path_sp=path_sp,
#                           psel=psel, fsel=fsel)
#
# @q.timer_verbose
# def compute_prop_wsrc_all(job_tag, traj, *,
#                           inv_type, gf, gt, wi, psel, fsel, eig):
#     inv_type_names = [ "light", "strange", ]
#     inv_type_name = inv_type_names[inv_type]
#     path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}"
#     path_sp = f"{job_tag}/psel-prop-wsrc-{inv_type_name}/traj-{traj}"
#     finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
#     sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", q.Coordinate([ 2, 2, 2, 4, ]))
#     for inv_acc in [ 2, 1, ]:
#         for p in wi:
#             idx, tslice, inv_type_p, inv_acc_p=p
#             if inv_type_p == inv_type and inv_acc_p == inv_acc:
#                 compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc,
#                         idx=idx, sfw=sfw, path_sp=path_sp,
#                         psel=psel, fsel=fsel, eig=eig,
#                         finished_tags=finished_tags)
#     sfw.close()
#     q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
#     # q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint ; wsnk.txt")))
#     q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
#     q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after=True)
#     # q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after=True)
#
# @q.timer
# def run_prop_wsrc(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt, get_psel, get_fsel, get_wi):
#     """
#     Can use `run_prop_wsrc_sparse` instead.
#     #
#     run_prop_wsrc_full(job_tag, traj, inv_type=0, get_gf=get_gf, get_eig=get_eig, get_gt=get_gt, get_wi=get_wi)
#     get_fsel, get_psel, get_fsel_prob, get_psel_prob = run_fsel_psel_from_wsrc_prop_full(job_tag, traj, get_wi=get_wi)
#     run_prop_wsrc_sparse(job_tag, traj, inv_type=0, get_gt=get_gt, get_psel=get_psel, get_fsel=get_fsel, get_wi=get_wi)
#     """
#     if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
#         return
#     if get_eig is None:
#         if inv_type == 0:
#             return
#         get_eig = lambda: None
#     inv_type_names = [ "light", "strange", ]
#     inv_type_name = inv_type_names[inv_type]
#     if get_load_path(f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}/geon-info.txt") is not None:
#         return
#     if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-{inv_type_name}"):
#         gf = get_gf()
#         gt = get_gt()
#         eig = get_eig()
#         psel = get_psel()
#         fsel = get_fsel()
#         assert fsel.is_containing(psel)
#         wi = get_wi()
#         compute_prop_wsrc_all(job_tag, traj,
#                               inv_type=inv_type, gf=gf, gt=gt, wi=wi,
#                               psel=psel, fsel=fsel, eig=eig)
#         q.release_lock()
