import numpy as np

from .c import *
from .utils import *
from .data import *
from .parallel import *

@timer
def build_corr_from_param_arr(
        param_arr,
        *,
        t_arr,
        n_ops,
        n_eigs=None,
        t_start_arr=None,
        t_size=None,
        atw_factor_arr=None,
        extra_state_sign_arr=None,
        extra_state_sign_t_start_arr=None,
        ):
    r"""
    param_arr = np.concatenate([ es.ravel(), cs.ravel() ], dtype=np.float64)
    es.shape == (n_eigs,)
    cs.shape == (n_eigs, n_ops,)
    #
    Also work if `jk_param_arr` is input in place of `param_arr`.
    #
    extra_state_sign_arr[ei] = extra_state_sign_arr[ei]
    state_sign_arr[ei] = (extra_state_sign_arr[ei]
                         * (es[ei] / abs(es[ei]))**((t_start_arr[ei] + extra_state_sign_t_start_arr) % 2)
                         )
    #
    corr_data.shape == (n_ops, n_ops, n_tslice,)
    corr_data[i, j, t] = (
        \sum_{ei} cs[ei, i] * cs[ei, j] * es[ei]**(t_arr[t] - t_start_arr[ei])
        * state_sign_arr[ei]
        )
    #
    if t_size is not None:
        corr_data[i, j, t] += (
            atw_factor_arr[i] * atw_factor_arr[j] *
            \sum_{ei} cs[ei, i] * cs[ei, j] * es[ei]**(t_size - t_arr[t] - t_start_arr[ei])
            * state_sign_arr[ei]
            )
    #
    return corr_data
    """
    if len(param_arr.shape) == 2:
        jk_param_arr = param_arr
        n_jk = jk_param_arr.shape[0]
        n_params = jk_param_arr.shape[1]
        if n_eigs is None:
            assert n_params % (n_ops + 1) == 0
            n_eigs = n_params // (n_ops + 1)
        else:
            assert n_params == n_eigs + n_eigs * n_ops
        n_tslice = t_arr.shape[0]
        assert jk_param_arr.shape == (n_jk, n_params,)
        assert t_arr.shape == (n_tslice,)
        jk_corr_data = np.zeros((n_jk, n_ops, n_ops, n_tslice,), dtype=np.float64)
        for jk_idx, param_arr in enumerate(jk_param_arr):
            jk_corr_data[jk_idx] = build_corr_from_param_arr(
                    param_arr, n_ops=n_ops, t_arr=t_arr,
                    t_start_arr=t_start_arr,
                    t_size=t_size,
                    atw_factor_arr=atw_factor_arr,
                    extra_state_sign_arr=extra_state_sign_arr,
                    extra_state_sign_t_start_arr=extra_state_sign_t_start_arr,
                    )
        return jk_corr_data
    assert len(t_arr.shape) == 1
    assert len(t_arr) >= 1
    assert len(param_arr.shape) == 1
    n_params = param_arr.shape[0]
    assert n_params % (n_ops + 1) == 0
    n_eigs = n_params // (n_ops + 1)
    param_arr = np.array(param_arr, dtype=np.float64)
    es_ini = param_arr[:n_eigs]
    es = np.zeros(n_eigs, dtype=np.float64)
    es[:] = es_ini
    cs = param_arr[n_eigs:].reshape(n_eigs, n_ops)
    t_start_arr_ini = t_start_arr
    t_start_arr = np.zeros(n_eigs, dtype=np.int32)
    if t_start_arr_ini is not None:
        t_start_arr[:] = t_start_arr_ini
    extra_state_sign_arr_ini = extra_state_sign_arr
    extra_state_sign_arr = np.ones(n_eigs, dtype=np.float64)
    if extra_state_sign_arr_ini is not None:
        extra_state_sign_arr[:] = extra_state_sign_arr_ini
    extra_state_sign_t_start_arr_ini = extra_state_sign_t_start_arr
    extra_state_sign_t_start_arr = np.zeros(n_eigs, dtype=np.int32)
    if extra_state_sign_t_start_arr_ini is not None:
        extra_state_sign_t_start_arr[:] = extra_state_sign_t_start_arr_ini
    state_sign_arr = extra_state_sign_arr * (es / abs(es))**((t_start_arr + extra_state_sign_t_start_arr) % 2)
    # corr_data[op1_idx, op2_idx, t_idx]
    corr_data = (
            cs[:, :, None, None] * cs[:, None, :, None]
            * es[:, None, None, None]**(t_arr - t_start_arr[:, None, None, None])
            * state_sign_arr[:, None, None, None]
            ).sum(0)
    if t_size is not None:
        atw_factor_arr_ini = atw_factor_arr
        atw_factor_arr = np.ones(n_ops, dtype=np.float64)
        if atw_factor_arr_ini is not None:
            atw_factor_arr[:] = atw_factor_arr_ini
        corr_data = corr_data + (
                cs[:, :, None, None] * cs[:, None, :, None]
                * atw_factor_arr[:, None, None]
                * atw_factor_arr[None, :, None]
                * es[:, None, None, None]**(t_size - t_arr - t_start_arr[:, None, None, None])
                * state_sign_arr[:, None, None, None]
                ).sum(0)
    return corr_data

@timer
def mk_data_set(
        *, n_jk=10, n_ops=4, n_eigs=4,
        t_arr=None,
        t_start_arr=None,
        extra_state_sign_arr=None,
        extra_state_sign_t_start_arr=None,
        t_size=None,
        atw_factor_arr=None,
        sigma=0.1,
        rng=None,
        ):
    r"""
    return param_arr, jk_corr_data, corr_data_sigma
    #
    es[ei] is float (eig with index ei)
    cs[ei, i] is float (coef with eig index ei and state index i)
    corr_data[i, j, t] = \sum_{ei} coefs[ei, i] * coefs[ei, j] * es[ei]**t
    jk_corr_data[jk, i, j, t] = corr_data[i, j, t] + corr_data_sigma[i, j, t] * N(0,1)
    #
    `t_arr`, `t_size`, `t_start_arr`, `atw_factor_arr`: see `build_corr_from_param_arr`
    """
    if rng is None:
        rng = RngState("mk_data_set-seed")
    #
    es = rng.u_rand_arr((n_eigs,)) * 1.2 - 0.3
    cs = rng.u_rand_arr((n_eigs, n_ops)) * 2.0 - 1.0
    param_arr = np.concatenate([ es.ravel(), cs.ravel() ], dtype=np.float64)
    #
    if t_arr is None:
        t_arr = np.arange(4)
    #
    n_tslice = t_arr.shape[0]
    assert t_arr.shape == (n_tslice,)
    #
    corr_data = build_corr_from_param_arr(
            param_arr, n_ops=n_ops, t_arr=t_arr,
            t_start_arr=t_start_arr,
            extra_state_sign_arr=extra_state_sign_arr,
            extra_state_sign_t_start_arr=extra_state_sign_t_start_arr,
            t_size=t_size,
            atw_factor_arr=atw_factor_arr,
            )
    shape = (n_ops, n_ops, n_tslice,)
    corr_data_sigma = np.zeros(shape, dtype=np.float64) + sigma
    #
    jk_corr_data = np.empty((n_jk,) + shape, dtype=np.float64)
    jk_corr_data = corr_data + rng.g_rand_arr((n_jk,) + shape) * corr_data_sigma
    #
    return param_arr, jk_corr_data, corr_data_sigma, t_arr

@timer
def sort_param_arr_free_eig(param_arr, n_ops, free_eig_idx_arr):
    """
    Adjust order of states with free eig parameter
    """
    if len(free_eig_idx_arr) == 0:
        return param_arr
    n_eigs = len(param_arr) // (n_ops + 1)
    assert len(param_arr) == n_eigs* (n_ops + 1)
    e_arr = param_arr[:n_eigs].copy()
    c_arr = param_arr[n_eigs:].reshape(n_eigs, n_ops).copy()
    free_e_arr = e_arr[free_eig_idx_arr]
    free_c_arr = c_arr[free_eig_idx_arr]
    free_e_arr, free_c_arr = zip(*sorted(zip(free_e_arr, free_c_arr), key=lambda p: -abs(p[0])))
    e_arr[free_eig_idx_arr] = free_e_arr
    c_arr[free_eig_idx_arr] = free_c_arr
    new_param_arr = np.concatenate([ e_arr, c_arr.ravel(), ], dtype=np.float64)
    return new_param_arr

@timer
def apply_eig_maximum(param_arr, eig_maximum_arr=None, free_eig_idx_arr=None):
    """
    return new param_arr (does not change original param_arr)
    param_arr has to be np.array
    should work for jk_param_arr as well
    #
    if eig_maximum_arr is None:
        return param_arr
    ...
    sel = abs(es) > eig_maximum_arr
    new_es[sel] = eig_maximum_arr[sel]**2 / es[sel]
    ...
    return new_param_arr
    """
    if eig_maximum_arr is None:
        return param_arr
    assert free_eig_idx_arr is not None
    new_param_arr = param_arr.copy()
    es = new_param_arr[..., free_eig_idx_arr]
    new_es = es.copy()
    max_es = es.copy()
    max_es[:] = eig_maximum_arr
    sel = abs(es) > max_es
    new_es[sel] = max_es[sel]**2 / es[sel]
    new_param_arr[..., free_eig_idx_arr] = new_es
    return new_param_arr

@timer
def mk_fcn(
        corr_data,
        corr_data_sigma,
        t_arr,
        t_start_arr,
        *,
        extra_state_sign_arr=None,
        extra_state_sign_t_start_arr=None,
        t_size=None,
        atw_factor_arr=None,
        eig_maximum_arr=None,
        free_eig_idx_arr=None,
        ):
    r"""
    Shape of inputs and parameters are like:
    `corr_data.shape == (n_ops, n_ops, t_len,)`
    `corr_data_sigma.shape == (n_ops, n_ops, t_len,)`
    `t_arr`, `t_size`, `t_start_arr`, `atw_factor_arr`: see `build_corr_from_param_arr`
    `eig_maximum_arr`, `free_eig_idx_arr`: see `apply_eig_maximum`
    `corr_data`: see `build_corr_from_param_arr`
    return fcn
    fcn(param_arr) => chisq, param_grad_arr
    `param_arr`: see `build_corr_from_param_arr`
    """
    import jax
    import jax.numpy as jnp
    if jnp.zeros(2, dtype=jnp.float64).dtype != jnp.float64:
        raise Exception(f"ERROR: double precision not available in JAX. Please set: 'export JAX_ENABLE_X64=True' to use double precision number here.")
    shape = corr_data.shape
    assert shape == corr_data_sigma.shape
    assert len(shape) == 3
    n_ops = shape[0]
    assert shape[1] == n_ops
    t_len = shape[2]
    t_arr = jnp.array(t_arr, dtype=jnp.int32)
    assert len(t_arr) == t_len
    t_start_arr = jnp.array(t_start_arr, dtype=jnp.int32)
    n_eigs = len(t_start_arr)
    extra_state_sign_arr_ini = extra_state_sign_arr
    extra_state_sign_arr = jnp.ones(n_eigs, dtype=jnp.float64)
    if extra_state_sign_arr_ini is not None:
        extra_state_sign_arr[:] = extra_state_sign_arr_ini
    extra_state_sign_t_start_arr_ini = extra_state_sign_t_start_arr
    extra_state_sign_t_start_arr = jnp.zeros(n_eigs, dtype=jnp.int32)
    if extra_state_sign_t_start_arr_ini is not None:
        extra_state_sign_t_start_arr[:] = extra_state_sign_t_start_arr_ini
    if atw_factor_arr is None:
        atw_factor_arr = jnp.ones(n_ops, dtype=jnp.float64)
    else:
        atw_factor_arr = jnp.array(atw_factor_arr, dtype=jnp.float64)
    if eig_maximum_arr is not None:
        assert free_eig_idx_arr is not None
        eig_maximum_arr = jnp.array(eig_maximum_arr, dtype=jnp.float64)
        free_eig_idx_arr = jnp.array(free_eig_idx_arr, dtype=jnp.int32)
        assert eig_maximum_arr.shape == free_eig_idx_arr.shape
    corr_avg = jnp.array(corr_data, dtype=jnp.float64)
    corr_sigma = jnp.array(corr_data_sigma, dtype=jnp.float64)
    def fcn_f(param_arr):
        assert len(param_arr.shape) == 1
        n_params = len(param_arr)
        assert n_params % (n_ops + 1) == 0
        n_eigs = n_params // (n_ops + 1)
        es = param_arr[:n_eigs]
        cs = param_arr[n_eigs:].reshape(n_eigs, n_ops)
        if eig_maximum_arr is not None:
            fes = es[free_eig_idx_arr]
            es = es.at[free_eig_idx_arr].set(
                    jnp.where(
                        abs(fes) > eig_maximum_arr,
                        eig_maximum_arr**2 / fes,
                        fes,
                        )
                    )
        state_sign_arr = extra_state_sign_arr * (es / abs(es))**((t_start_arr + extra_state_sign_t_start_arr) % 2)
        # corr[op1_idx, op2_idx, t_idx]
        corr = (cs[:, :, None, None] * cs[:, None, :, None]
                * es[:, None, None, None]**(t_arr - t_start_arr[:, None, None, None])
                * state_sign_arr[:, None, None, None]
                ).sum(0)
        if t_size is not None:
            corr = corr + (
                    cs[:, :, None, None] * cs[:, None, :, None]
                    * atw_factor_arr[:, None, None]
                    * atw_factor_arr[None, :, None]
                    * es[:, None, None, None]**(t_size - t_arr - t_start_arr[:, None, None, None])
                    * state_sign_arr[:, None, None, None]
                    ).sum(0)
        corr = (corr - corr_avg) / corr_sigma
        chisqs = corr * corr
        chisq = chisqs.sum()
        return chisq
    fcn_f_jit = jax.jit(fcn_f)
    fcn_fg_jit = jax.jit(jax.value_and_grad(fcn_f))
    @timer
    def fcn(param_arr, requires_grad=True):
        param_arr = jnp.array(param_arr, dtype=jnp.float64)
        if requires_grad:
            chisq, grad = fcn_fg_jit(param_arr)
            return chisq.item(), np.array(grad, dtype=np.float64)
        else:
            chisq = fcn_f_jit(param_arr)
            return chisq.item()
    return fcn

### -------------------

@timer
def minimize(fcn, n_step=10, step_size=1e-2, *, param_arr):
    fname = get_fname()
    chisq_pre = None
    for i in range(n_step):
        chisq, param_grad_arr = fcn(param_arr)
        grad_norm = np.linalg.norm(param_grad_arr)
        displayln_info(1, f"{fname}: step={i} chisq={chisq} grad_norm={grad_norm} param_arr={param_arr}")
        if chisq_pre is not None:
            if chisq > chisq_pre:
                displayln_info(0, f"{fname}: early stop step={i} step_size={step_size} chisq={chisq_pre} grad_norm={grad_norm}")
                return param_arr_pre, i
        param_arr_pre = param_arr
        chisq_pre = chisq
        param_arr = param_arr - step_size / grad_norm * param_grad_arr
    displayln_info(0, f"{fname}: step={n_step} step_size={step_size} chisq={chisq_pre} grad_norm={grad_norm}")
    return param_arr_pre, n_step

@timer
def adaptive_minimize(fcn, step_size_list, n_step=10, max_total_steps=10000, *, param_arr):
    fname = get_fname()
    idx = 0
    total_steps = 0
    total_steps_pre = 0
    while True:
        step_size = step_size_list[idx]
        param_arr, i_step = minimize(fcn, n_step=n_step, step_size=step_size, param_arr=param_arr)
        total_steps += i_step
        if total_steps - total_steps_pre > 1000:
            displayln_info(0, f"{fname}: {step_size:.8f} total_steps={total_steps}")
            total_steps_pre = total_steps
        if total_steps > max_total_steps:
            return param_arr
        if i_step == n_step:
            idx = max(idx - 1, 0)
        else:
            idx = idx + 1
            if idx == len(step_size_list):
                return param_arr

@timer
def minimize_scipy(fcn, *, param_arr, fixed_param_mask=None, minimize_kwargs=None):
    import scipy
    fname = get_fname()
    n_params = len(param_arr)
    if fixed_param_mask is None:
        fixed_param_mask = np.zeros(n_params, dtype=bool)
    p_fixed = param_arr[fixed_param_mask]
    free_param_mask = ~fixed_param_mask
    @timer
    def c_fcn(p_free):
        p_all = np.empty(n_params, dtype=np.float64)
        p_all[fixed_param_mask] = p_fixed
        p_all[free_param_mask] = p_free
        chisq = fcn(p_all, requires_grad=False)
        return chisq
    @timer
    def c_fcn_grad(p_free):
        p_all = np.empty(n_params, dtype=np.float64)
        p_all[fixed_param_mask] = p_fixed
        p_all[free_param_mask] = p_free
        chisq, grad_all = fcn(p_all, requires_grad=True)
        return grad_all[free_param_mask]
    p_free = param_arr[free_param_mask]
    fcn_initial = c_fcn(p_free)
    minimize_kwargs_default = dict(
            options=dict(maxiter=1e3),
            method="BFGS",
            )
    if minimize_kwargs is None:
        minimize_kwargs = minimize_kwargs_default
    else:
        assert isinstance(minimize_kwargs, dict)
        minimize_kwargs = minimize_kwargs_default | minimize_kwargs
    res = scipy.optimize.minimize(
            c_fcn,
            p_free,
            jac=c_fcn_grad,
            **minimize_kwargs,
            )
    p_free_mini = res.x
    fcn_final = c_fcn(p_free_mini)
    displayln_info(0, f"{fname}: fun={res.fun} ; grad_norm={np.linalg.norm(res.jac)}")
    displayln_info(0, f"{fname}: success={res.success} ; message={res.message} ; nfev={res.nfev} ; njev={res.njev}")
    if fcn_final <= fcn_initial:
        param_arr_mini = param_arr.copy()
        param_arr_mini[free_param_mask] = p_free_mini
        return param_arr_mini
    else:
        displayln_info(0, f"{fname}: return initial parameter instead due to: fcn_initial={fcn_initial} < fcn_final={fcn_final} .")
        return param_arr

### -----------------

def mp_initializer():
    import qlat_utils as q
    q.set_verbose_level(-1)
    import jax
    def f(x):
        return x*x
    gf = jax.jit(jax.value_and_grad(f))
    v1, v2, = gf(1.0)
    assert (v1.item(), v2.item()) == (1.0, 2.0,)

def jk_mini_task_in_fit_eig_coef(kwargs):
    fname = get_fname()
    def f(*,
          jk_idx,
          corr_data,
          corr_data_err,
          t_arr,
          t_start_arr,
          extra_state_sign_arr,
          extra_state_sign_t_start_arr,
          t_size,
          atw_factor_arr,
          eig_maximum_arr,
          free_eig_idx_arr,
          param_arr_mini,
          n_step_mini_jk,
          r_amp,
          all_eig_mask,
          fixed_eig_mask,
          free_eig_mask,
          fixed_coef_eig_mask,
          minimize_kwargs,
          is_sorting_eig_state,
          rng_seed,
          verbose_level,
          ):
        set_verbose_level(verbose_level)
        rng = RngState(rng_seed)
        n_params = len(param_arr_mini)
        n_ops = corr_data.shape[0]
        fcn = mk_fcn(corr_data, corr_data_err, t_arr, t_start_arr,
                     extra_state_sign_arr=extra_state_sign_arr,
                     extra_state_sign_t_start_arr=extra_state_sign_t_start_arr,
                     t_size=t_size,
                     atw_factor_arr=atw_factor_arr,
                     eig_maximum_arr=eig_maximum_arr,
                     free_eig_idx_arr=free_eig_idx_arr,
                     )
        rand_update_mask = (~all_eig_mask) & (~fixed_coef_eig_mask)
        def display_param_arr(param_arr, mask=None, verbose_level=0):
            fcn_v, grad = fcn(param_arr)
            grad_norm = np.linalg.norm(grad)
            if mask is not None:
                grad_masked = grad[~mask]
            else:
                grad_masked = grad
            grad_norm_masked = np.linalg.norm(grad_masked)
            displayln_info(verbose_level, f"{fname}: fcn={fcn_v:.5E} grad_norm={grad_norm:.5E} grad_norm_masked={grad_norm_masked:.5E}")
            eigs = param_arr[all_eig_mask]
            grad_eigs = grad[all_eig_mask]
            grad_eigs_norm = np.linalg.norm(grad_eigs)
            eg_arr = np.stack([ eigs, grad_eigs, ]).T
            important_eg_arr = eg_arr[abs(grad_eigs) > grad_eigs_norm / 10]
            displayln_info(verbose_level, f"{fname}: eigs and grad arr=\n{important_eg_arr}")
        def rand_update(param_arr):
            param_arr = param_arr.copy()
            param_arr[rand_update_mask] = (
                    param_arr[rand_update_mask]
                    + (rng.u_rand_arr(n_params)[rand_update_mask] - 0.5) * r_amp)
            return param_arr
        param_arr = param_arr_mini.copy()
        display_param_arr(param_arr, mask=fixed_eig_mask, verbose_level=0)
        displayln_info(0, f"{fname}: mini fcn (fixed all eigs)")
        for i in range(n_step_mini_jk):
            param_arr = rand_update(param_arr)
            param_arr = minimize_scipy(fcn, param_arr=param_arr,
                                       fixed_param_mask=all_eig_mask | fixed_coef_eig_mask,
                                       minimize_kwargs=minimize_kwargs)
            vl = 1
            if i == n_step_mini_jk - 1:
                vl = 0
            display_param_arr(param_arr, mask=all_eig_mask, verbose_level=vl)
        displayln_info(0, f"{fname}: initial free_eig_arr={param_arr[free_eig_mask].tolist()}")
        displayln_info(0, f"{fname}: mini fcn (free eigs selected by free_eig_idx_arr)")
        for i in range(n_step_mini_jk):
            param_arr = rand_update(param_arr)
            param_arr = minimize_scipy(fcn, param_arr=param_arr,
                                       fixed_param_mask=fixed_eig_mask | fixed_coef_eig_mask,
                                       minimize_kwargs=minimize_kwargs)
            if is_sorting_eig_state:
                param_arr = sort_param_arr_free_eig(param_arr, n_ops, free_eig_idx_arr)
            param_arr = apply_eig_maximum(param_arr, eig_maximum_arr, free_eig_idx_arr)
            if is_sorting_eig_state:
                # Assuming `eig_maximum_arr` is sorted.
                param_arr = sort_param_arr_free_eig(param_arr, n_ops, free_eig_idx_arr)
            displayln_info(0, f"{fname}: iter={i} free_eig_arr={param_arr[free_eig_mask].tolist()}")
            vl = 1
            if i == n_step_mini_jk - 1:
                vl = 0
            display_param_arr(param_arr, mask=fixed_eig_mask, verbose_level=vl)
        chisq, param_grad = fcn(param_arr)
        set_verbose_level(-1)
        return chisq, param_arr, param_grad
    return f(**kwargs)

@timer_verbose
def mk_mp_pool(n_proc=None):
    """
    return `mp_pool`
    #
    Usage of `mp_pool`
    mp_map = mp_pool.imap
    """
    if n_proc is None:
        n_proc = get_q_num_mp_processes()
    assert isinstance(n_proc, int)
    mp_pool_n_proc = n_proc
    import multiprocessing
    mp_pool = multiprocessing.get_context('spawn').Pool(mp_pool_n_proc, initializer=mp_initializer)
    mp_map = mp_pool.imap
    assert list(mp_map(np.sin, range(n_proc))) == list(map(np.sin, range(n_proc)))
    return mp_pool

@timer_verbose
def close_mp_pool(mp_pool):
    mp_pool.close()

@timer_verbose
def fit_eig_coef(jk_corr_data,
                 *,
                 t_arr,
                 e_arr,
                 t_start_arr=None,
                 extra_state_sign_arr=None,
                 extra_state_sign_t_start_arr=None,
                 t_size=None,
                 atw_factor_arr=None,
                 eig_maximum_arr=None,
                 c_arr=None,
                 op_norm_fac_arr=None,
                 free_eig_idx_arr=None,
                 fixed_coef_eig_idx_arr=None,
                 n_step_mini_avg=10,
                 n_step_mini_jk=5,
                 minimize_kwargs=None,
                 r_amp=1e-3,
                 diag_err_scale_factor=1.0,
                 off_diag_err_scale_factor=1.0,
                 rng_seed_list=None,
                 mp_pool=None,
                 ):
    r"""
    return res
    #
    res['jk_chisq'] = jk_chisq
    res['jk_param_arr'] = jk_param_arr
    res['jk_param_grad_arr'] = jk_param_grad_arr # chisq grad with respect to param_arr
    res['jk_param_for_scaled_corr_arr'] = jk_param_for_scaled_corr_arr
    res['jk_param_grad_for_scaled_corr_arr'] = jk_param_grad_for_scaled_corr_arr
    #
    C^{(jk_idx)}_{i,j}(t_arr[t]) == jk_corr_data[jk_idx, i, j, t]
    #
    jk_corr_data.shape == (n_jk, n_ops, n_ops, n_tslice,)
    jk_corr_data.dtype == np.float64
    #
    free_eig_arr = e_arr[free_eig_idx_arr]
    # by default, all eigs are fixed
    #
    e_arr.shape == (n_eigs,)
    c_arr.shape == (n_eigs, n_ops,)
    #
    C^{(jk_idx)}_{i,j}(t_arr[t]) ~ \sum_{n} c_arr[n, i] * c_arr[n, j] * e_arr[n]**(t_arr[t] - t_start_arr[n]) * state_sign_arr[n]
    state_sign_arr = extra_state_sign_arr * (es / abs(es))**((t_start_arr + extra_state_sign_t_start_arr) % 2)
    #
    jk_param_arr_mini.shape == (n_jk, n_params)
    param_arr == np.concatenate([ e_arr, c_arr.ravel(), ], dtype=np.float64)
    #
    diag_err_scale_factor should be 1.0
    off_diag_err_scale_factor should be np.sqrt(2) if jk_corr_data has been symmetrized.
    #
    `eig_maximum_arr` should be of same shape as `free_eig_idx_arr` will constrain all the free eigs to be smaller than this eig value (None means no constraint)
    #
    rng_seed_list=[ f"fit-eig-coef-seed-{i}" for i in range(32) ]
    mp_pool = mk_mp_pool(n_proc)
    """
    fname = get_fname()
    #
    verbose_level = get_verbose_level()
    #
    assert len(jk_corr_data.shape) == 4
    t_arr = np.array(t_arr, dtype=np.int32)
    e_arr = np.array(e_arr, dtype=np.float64)
    assert len(t_arr.shape) == 1
    assert len(e_arr.shape) == 1
    #
    jk_corr_data = jk_corr_data.copy()
    #
    n_jk = jk_corr_data.shape[0]
    n_ops = jk_corr_data.shape[1]
    n_tslice = t_arr.shape[0]
    n_eigs = e_arr.shape[0]
    #
    assert jk_corr_data.shape == (n_jk, n_ops, n_ops, n_tslice,)
    assert t_arr.shape == (n_tslice,)
    assert e_arr.shape == (n_eigs,)
    #
    if free_eig_idx_arr is None:
        free_eig_idx_arr = np.array([], dtype=np.int32)
    else:
        free_eig_idx_arr = np.array(free_eig_idx_arr, dtype=np.int32)
    #
    if eig_maximum_arr is not None:
        eig_maximum_arr = np.array(eig_maximum_arr, dtype=np.float64)
        assert eig_maximum_arr.shape == free_eig_idx_arr.shape
    #
    if fixed_coef_eig_idx_arr is None:
        fixed_coef_eig_idx_arr = np.array([], dtype=np.int32)
    else:
        fixed_coef_eig_idx_arr = np.array(fixed_coef_eig_idx_arr, dtype=np.int32)
    #
    n_free_eigs = free_eig_idx_arr.shape[0]
    n_fixed_eigs = n_eigs - n_free_eigs
    n_fixed_coef_eigs = fixed_coef_eig_idx_arr.shape[0]
    #
    assert free_eig_idx_arr.shape == (n_free_eigs,)
    assert fixed_coef_eig_idx_arr.shape == (n_fixed_coef_eigs,)
    #
    if c_arr is None:
        c_arr = np.zeros((n_eigs, n_ops,), dtype=np.float64)
    else:
        assert c_arr.shape == (n_eigs, n_ops,)
    #
    t_start_arr_ini = t_start_arr
    t_start_arr = np.zeros(n_eigs, dtype=np.int32)
    if t_start_arr_ini is not None:
        t_start_arr[:] = t_start_arr_ini
    #
    atw_factor_arr_ini = atw_factor_arr
    atw_factor_arr = np.ones(n_ops, dtype=np.float64)
    if atw_factor_arr_ini is not None:
        atw_factor_arr[:] = atw_factor_arr_ini
    #
    op_idx_arr = np.arange(n_ops)
    op_norm_fac_arr_ini = op_norm_fac_arr
    op_norm_fac_arr = 1 / np.sqrt(abs(jk_corr_data[0, op_idx_arr, op_idx_arr, 0]))
    if op_norm_fac_arr_ini is not None:
        op_norm_fac_arr[:] = op_norm_fac_arr_ini
    #
    jk_corr_data = op_norm_fac_arr[:, None, None] * op_norm_fac_arr[None, :, None] * jk_corr_data
    #
    if rng_seed_list is None:
        rng_seed_list = [ "fit_eig_coef-seed-param", ]
    elif isinstance(rng_seed_list, str):
        rng_seed_list = [ rng_seed_list, ]
    else:
        assert len(rng_seed_list) >= 1
        for v in rng_seed_list:
            assert isinstance(v, str)
    #
    if n_step_mini_jk == 0:
        mp_pool = None
    #
    is_close_pool = False
    if mp_pool is None:
        mp_map = map
        mp_pool_n_proc = 1
    elif isinstance(mp_pool, int):
        mp_pool_n_proc = mp_pool
        mp_pool = mk_mp_pool(mp_pool_n_proc)
        is_close_pool = True
        mp_map = mp_pool.imap
    else:
        mp_map = mp_pool.imap
    #
    corr_data, corr_data_err = g_jk_avg_err(jk_corr_data)
    #
    is_finite_sel = np.isfinite(corr_data)
    jk_corr_data[:, ~is_finite_sel] = 0.0
    corr_data[~is_finite_sel] = 0.0
    corr_data_err[~is_finite_sel] = np.inf
    #
    is_zero_err_sel = corr_data_err == 0.0
    corr_data_err[is_zero_err_sel] = 1.0
    #
    # corr_data_err[op1_idx, op2_idx, t_idx]
    op_idx_arr = np.arange(n_ops)
    op_op_diag_sel = op_idx_arr[:, None] == op_idx_arr[None, :]
    corr_data_err[op_op_diag_sel] *= diag_err_scale_factor
    corr_data_err[~op_op_diag_sel] *= off_diag_err_scale_factor
    #
    c_arr = c_arr * op_norm_fac_arr
    #
    param_arr_initial = np.concatenate([ e_arr, c_arr.ravel(), ], dtype=np.float64)
    chisq_initial = np.inf
    #
    n_params = len(param_arr_initial)
    #
    all_eig_mask = np.arange(n_params) < n_eigs
    fixed_eig_mask = all_eig_mask.copy()
    fixed_eig_mask[free_eig_idx_arr] = False
    #
    c_mask = np.full((n_eigs, n_ops,), False, dtype=bool)
    c_mask[fixed_coef_eig_idx_arr] = True
    fixed_coef_eig_mask = np.full(n_params, False, dtype=bool)
    fixed_coef_eig_mask[n_eigs:] = c_mask.ravel()
    free_eig_mask = (~fixed_eig_mask) & all_eig_mask
    #
    param_arr_mini = param_arr_initial.copy()
    chisq_mini = chisq_initial
    #
    def mk_kwargs(jk_idx, is_sorting_eig_state, rng_seed, verbose_level):
        kwargs = dict(
                jk_idx=jk_idx,
                corr_data=jk_corr_data[jk_idx],
                corr_data_err=corr_data_err,
                t_arr=t_arr,
                t_start_arr=t_start_arr,
                extra_state_sign_arr=extra_state_sign_arr,
                extra_state_sign_t_start_arr=extra_state_sign_t_start_arr,
                t_size=t_size,
                atw_factor_arr=atw_factor_arr,
                eig_maximum_arr=eig_maximum_arr,
                free_eig_idx_arr=free_eig_idx_arr,
                param_arr_mini=param_arr_mini,
                n_step_mini_jk=n_step_mini_jk,
                r_amp=r_amp,
                all_eig_mask=all_eig_mask,
                fixed_eig_mask=fixed_eig_mask,
                free_eig_mask=free_eig_mask,
                fixed_coef_eig_mask=fixed_coef_eig_mask,
                minimize_kwargs=minimize_kwargs,
                is_sorting_eig_state=is_sorting_eig_state,
                rng_seed=rng_seed,
                verbose_level=verbose_level,
                )
        return kwargs
    #
    displayln_info(0, f"{fname}: initial free_eig_arr={param_arr_mini[free_eig_mask].tolist()}")
    displayln_info(0, f"{fname}: mini avg with all rng_seed_list")
    v_list = []
    for idx, v in enumerate(mp_map(jk_mini_task_in_fit_eig_coef,
                                   [ mk_kwargs(
                                       jk_idx=0,
                                       is_sorting_eig_state=True,
                                       rng_seed=rng_seed,
                                       verbose_level=verbose_level if idx == 0 else -1,
                                       )
                                    for idx, rng_seed in enumerate(rng_seed_list)
                                    ])):
        set_verbose_level(verbose_level)
        v_list.append(v)
        chisq, param_arr, param_grad_arr, = v
        displayln_info(0, f"{fname}: map: rs_idx={idx} ; chisq={chisq} ; free_eig_arr={param_arr[free_eig_mask].tolist()} ; rng_seed='{rng_seed_list[idx]}'")
        if eig_maximum_arr is not None:
            displayln_info(0, f"{fname}: map: rs_idx={idx} ; free_eig_arr/eig_maximum_arr={(param_arr[free_eig_idx_arr]/eig_maximum_arr).tolist()}")
    #
    rng_seed_mini = rng_seed_list[0]
    for idx, v in enumerate(v_list):
        chisq, param_arr, param_grad_arr, = v
        if chisq < chisq_mini:
            chisq_mini = chisq
            rng_seed_mini = rng_seed_list[idx]
            param_arr_mini = param_arr
    #
    displayln_info(0, f"{fname}: chisq_mini={chisq_mini} ; rng_seed_mini='{rng_seed_mini}'")
    displayln_info(0, f"{fname}: avg mini free_eig_arr={param_arr_mini[free_eig_mask].tolist()}")
    if eig_maximum_arr is not None:
        displayln_info(0, f"{fname}: map: rs_idx={idx} ; free_eig_arr/eig_maximum_arr={(param_arr_mini[free_eig_idx_arr]/eig_maximum_arr).tolist()}")
    #
    displayln_info(0, f"{fname}: mini all jk samples")
    jk_chisq = []
    jk_param_arr = []
    jk_param_grad_arr = []
    for idx, v in enumerate(mp_map(jk_mini_task_in_fit_eig_coef,
                                   [ mk_kwargs(
                                       jk_idx=jk_idx,
                                       is_sorting_eig_state=False,
                                       rng_seed=f"{rng_seed_list[0]}/jk_mini_task/{jk_idx}",
                                       verbose_level=verbose_level if jk_idx == 0 else -1,
                                       )
                                    for jk_idx in range(n_jk)
                                    ])):
        set_verbose_level(verbose_level)
        chisq, param_arr, param_grad_arr, = v
        jk_chisq.append(chisq)
        jk_param_grad_arr.append(param_grad_arr)
        jk_param_arr.append(param_arr)
        if n_step_mini_jk != 0:
            displayln_info(0, f"{fname}: map: jk_idx={idx} ; chisq={chisq} ; free_eig_arr={param_arr[free_eig_mask].tolist()}")
            if eig_maximum_arr is not None:
                displayln_info(0, f"{fname}: map: jk_idx={idx} ; free_eig_arr/eig_maximum_arr={(param_arr[free_eig_idx_arr]/eig_maximum_arr).tolist()}")
    if is_close_pool:
        mp_pool.close()
    #
    jk_chisq = np.array(jk_chisq, dtype=np.float64)
    jk_param_arr = np.array(jk_param_arr, dtype=np.float64)
    jk_param_grad_arr = np.array(jk_param_grad_arr, dtype=np.float64)
    #
    jk_param_grad_for_scaled_corr_arr = jk_param_grad_arr.copy()
    jk_param_for_scaled_corr_arr = jk_param_arr.copy()
    #
    jk_e_arr = jk_param_arr[:, :n_eigs].copy()
    jk_c_arr = jk_param_arr[:, n_eigs:].reshape(n_jk, n_eigs, n_ops,).copy()
    jk_c_arr = jk_c_arr / op_norm_fac_arr
    jk_param_arr[:, :n_eigs] = jk_e_arr
    jk_param_arr[:, n_eigs:] = jk_c_arr.reshape(n_jk, n_eigs * n_ops)
    #
    jk_e_grad_arr = jk_param_grad_arr[:, :n_eigs].copy()
    jk_c_grad_arr = jk_param_grad_arr[:, n_eigs:].reshape(n_jk, n_eigs, n_ops,).copy()
    jk_c_grad_arr = jk_c_grad_arr / op_norm_fac_arr
    jk_param_grad_arr[:, :n_eigs] = jk_e_grad_arr
    jk_param_grad_arr[:, n_eigs:] = jk_c_grad_arr.reshape(n_jk, n_eigs * n_ops)
    #
    res = dict()
    res['jk_chisq'] = jk_chisq
    res['jk_param_arr'] = jk_param_arr
    res['jk_param_grad_arr'] = jk_param_grad_arr
    res['jk_param_for_scaled_corr_arr'] = jk_param_for_scaled_corr_arr
    res['jk_param_grad_for_scaled_corr_arr'] = jk_param_grad_for_scaled_corr_arr
    res['jk_e_arr'] = jk_e_arr
    res['jk_c_arr'] = jk_c_arr
    res['jk_e_grad_arr'] = jk_e_grad_arr
    res['jk_c_grad_arr'] = jk_c_grad_arr
    res['options'] = dict(
            t_arr=t_arr,
            n_ops=n_ops,
            n_eigs=n_eigs,
            t_start_arr=t_start_arr,
            t_size=t_size,
            atw_factor_arr=atw_factor_arr,
            extra_state_sign_arr=extra_state_sign_arr,
            extra_state_sign_t_start_arr=extra_state_sign_t_start_arr,
            )
    #
    displayln_info(0, f"{fname} finished")
    return res
