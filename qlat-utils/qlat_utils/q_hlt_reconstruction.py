__all__ = [
        'mk_hlt_params',
        #
        'delta_from_g',
        #
        'aa_from_g_via_sum',
        'normalization_constraint_via_sum',
        'ww_from_g_via_sum',
        'ww_from_g_wgrad_via_sum',
        'mk_g_t_arr_via_sum',
        #
        'aa_from_g',
        'normalization_constraint',
        'ww_from_g',
        'ww_from_g_wgrad',
        'mk_g_t_arr',
        #
        'get_f_e_weight_log',
        ]

import jax
import jax.numpy as jnp

import numpy as np
import scipy.integrate as integrate

import qlat_utils as q

if jnp.zeros(2, dtype=jnp.float64).dtype != jnp.float64:
    raise Exception(f"ERROR: double precision not available in JAX. Please set: 'export JAX_ENABLE_X64=True' to use double precision number here.")

##############

def mk_hlt_params():
    """
    t_arr.shape == (t_size,)
    e_arr.shape == (n_energies,)
    type(cov.shape) == float or cov.shape == (t_size,) or cov.shape == (t_size, t_size,)
    """
    params = {}
    # mandatory params
    params["f_delta_target"] = None
    params["t_arr"] = None
    params["cov"] = None
    # mandatory params for via sum functions
    params["e_arr"] = None
    # optional params
    params["e0"] = 0.0
    params["ee_max"] = np.inf # only used for functions use integration with e
    params["lambda"] = 1.0
    params["alpha"] = -0.01 # override by "f_e_weight_log", see `get_f_e_weight_log`
    params["f_e_weight_log"] = None
    params["tt_size"] = None
    params["atw_factor"] = 1.0 # only has effects if "tt_size" is not None
    params["minimization_iter_max"] = 1
    params["g_t_arr_init"] = None
    params["does_have_constraint"] = False
    return params

def get_f_e_weight_log(params):
    f = params["f_e_weight_log"]
    if f is not None:
        return f
    alpha = params["alpha"]
    def f(e):
        return alpha * e
    return f

##############

@jax.jit
def delta_from_g(g_t_arr, t_arr, e_arr):
    """
    return delta
    g_t_arr.shape == (..., t_size,)
    delta.shape[-1] == (..., n_energies,)
    """
    delta = (jnp.exp(-e_arr[:, None] * t_arr) * g_t_arr[..., None, :]).sum(-1)
    return delta

############## with summation ##############

def aa_from_g_via_sum(g_t_arr, params):
    """
    g_t_arr.shape == (t_size,)
    delta_target.shape == (n_energies,)
    """
    t_arr = params["t_arr"]
    e_arr = params["e_arr"]
    f_delta_target = params["f_delta_target"]
    tt_size = params["tt_size"]
    atw_factor = params["atw_factor"]
    delta_target = f_delta_target(e_arr)
    f_e_weight_log = get_f_e_weight_log(params)
    e0 = params["e0"]
    e_w = np.exp(f_e_weight_log(e_arr))
    e_w[e_arr < e0] = 0
    e_w[0] *= (e_arr[1] - e_arr[0]) / 2
    e_w[-1] *= (e_arr[-1] - e_arr[-2]) / 2
    e_w[1:-1] *= (e_arr[2:] - e_arr[:-2]) / 2
    delta = delta_from_g(g_t_arr, t_arr, e_arr)
    if tt_size is not None:
        delta = delta + atw_factor * delta_from_g(g_t_arr, tt_size - t_arr, e_arr)
    f_e = (delta - delta_target)**2 * e_w
    return f_e.sum()

def get_cov_term(g_t_arr, cov):
    t_size = len(g_t_arr)
    assert g_t_arr.shape == (t_size,)
    if len(cov.shape) < 2:
        if len(cov.shape) == 1:
            assert cov.shape == (t_size,)
        ee = (cov * g_t_arr * g_t_arr).sum()
    else:
        assert cov.shape == (t_size, t_size,)
        ee = (cov * g_t_arr[:, None] * g_t_arr).sum()
    return ee

def normalization_constraint_via_sum(g_t_arr, params):
    """
    return new_g_t_arr, constraint_penalty
    return a new g_t_arr which satisfy constraint and penalty reflect the deviation
    """
    does_have_constraint = params["does_have_constraint"]
    if not does_have_constraint:
        return g_t_arr, 0.0
    t_arr = params["t_arr"]
    e_arr = params["e_arr"]
    f_delta_target = params["f_delta_target"]
    tt_size = params["tt_size"]
    atw_factor = params["atw_factor"]
    delta_target = f_delta_target(e_arr)
    f_e_weight_log = get_f_e_weight_log(params)
    e0 = params["e0"]
    e_w = np.ones_like(e_arr)
    e_w[e_arr < e0] = 0
    e_w[0] *= (e_arr[1] - e_arr[0]) / 2
    e_w[-1] *= (e_arr[-1] - e_arr[-2]) / 2
    e_w[1:-1] *= (e_arr[2:] - e_arr[:-2]) / 2
    g1_t_arr = jnp.ones_like(g_t_arr)
    delta = delta_from_g(g_t_arr, t_arr, e_arr)
    delta1 = delta_from_g(g1_t_arr, t_arr, e_arr)
    if tt_size is not None:
        delta = delta + atw_factor * delta_from_g(g_t_arr, tt_size - t_arr, e_arr)
        delta1 = delta1 + atw_factor * delta_from_g(g1_t_arr, tt_size - t_arr, e_arr)
    s_target = (delta_target * e_w).sum()
    s = (delta * e_w).sum()
    s1 = (delta1 * e_w).sum()
    coef = (s_target - s) / s1
    new_g_t_arr = g_t_arr + coef * g1_t_arr
    constraint_penalty = 10.0 * coef**2
    return new_g_t_arr, constraint_penalty

def ww_from_g_via_sum(g_t_arr, params):
    """
    g_t_arr.shape == (t_size,)
    delta_target.shape == (n_energies,)
    """
    cov = params["cov"]
    t_arr = params["t_arr"]
    new_g_t_arr, constraint_penalty = normalization_constraint_via_sum(g_t_arr, params)
    g0_t_arr = jnp.zeros(new_g_t_arr.shape, jnp.float64)
    aa = aa_from_g_via_sum(new_g_t_arr, params) / aa_from_g_via_sum(g0_t_arr, params)
    ee = get_cov_term(new_g_t_arr, cov) * params["lambda"]
    return aa + ee + constraint_penalty

ww_from_g_wgrad_via_sum = jax.value_and_grad(ww_from_g_via_sum)

def mk_g_t_arr_optimization_fcn_via_sum(params):
    def fcn(g_t_arr, requires_grad=True):
        if requires_grad:
            return ww_from_g_wgrad_via_sum(g_t_arr, params)
        else:
            return ww_from_g_via_sum(g_t_arr, params)
    return fcn

@q.timer
def mk_g_t_arr_via_sum(params):
    fcn = mk_g_t_arr_optimization_fcn_via_sum(params)
    t_arr = params["t_arr"]
    g_t_arr = params["g_t_arr_init"]
    if g_t_arr is None:
        g_t_arr = np.zeros_like(t_arr, dtype=np.float64)
    for i in range(params["minimization_iter_max"]):
        g_t_arr = q.q_fit_corr.minimize_scipy(fcn, param_arr=g_t_arr)
    return g_t_arr

############## with integration ##############

def build_hlt_aa_mat(params):
    tag = "aa_mat"
    if tag in params:
        return params[tag]
    f_e_weight_log = get_f_e_weight_log(params)
    e0 = params["e0"]
    t_arr = params["t_arr"]
    tt_size = params["tt_size"]
    atw_factor = params["atw_factor"]
    ee_max = params["ee_max"]
    t_size = len(t_arr)
    aa_mat = np.zeros((t_size, t_size,), dtype=np.float64)
    aa_values = {}
    def compute_t_sum(t_sum):
        t_sum = int(t_sum)
        if t_sum not in aa_values:
            def f(e):
                return np.exp(f_e_weight_log(e) - e * t_sum)
            v = integrate.quad(f, e0, ee_max)[0]
            aa_values[t_sum] = v
        return aa_values[t_sum]
    for t1_idx, t1 in enumerate(t_arr):
        for t2_idx, t2 in enumerate(t_arr):
            v = compute_t_sum(t1 + t2)
            if tt_size is not None:
                v = v + atw_factor * compute_t_sum(t1 + tt_size - t2)
                v = v + atw_factor * compute_t_sum(tt_size - t1 + t2)
                v = v + atw_factor * atw_factor * compute_t_sum(tt_size - t1 + tt_size - t2)
            aa_mat[t1_idx, t2_idx] = v
    aa_mat = jnp.array(aa_mat)
    params[tag] = aa_mat
    return aa_mat

def build_hlt_f_vec(params):
    tag = "f_vec"
    if tag in params:
        return params[tag]
    f_e_weight_log = get_f_e_weight_log(params)
    e0 = params["e0"]
    t_arr = params["t_arr"]
    tt_size = params["tt_size"]
    atw_factor = params["atw_factor"]
    f_delta_target = params["f_delta_target"]
    ee_max = params["ee_max"]
    t_size = len(t_arr)
    f_vec = np.zeros(t_size, dtype=np.float64)
    f_values = {}
    def compute_t(t):
        t = int(t)
        if t not in f_values:
            def f(e):
                return -2 * f_delta_target(e) * np.exp(f_e_weight_log(e) - e * t)
            v = integrate.quad(f, e0, ee_max)[0]
            f_values[t] = v
        return f_values[t]
    for t_idx, t in enumerate(t_arr):
        v = compute_t(t)
        if tt_size is not None:
            v = v + atw_factor * compute_t(tt_size - t)
        f_vec[t_idx] = v
    f_vec = jnp.array(f_vec)
    params[tag] = f_vec
    return f_vec

def build_hlt_fc_vec(params):
    """
    return norm, fc_vec
    where
    norm = integrate.quad(f_delta_target, e0, ee_max)[0]
    normalization constraint can be expressed as:
    norm == (fc_vec * g_t_arr).sum()
    """
    tag = "fc_vec"
    if tag in params:
        return params[tag]
    e0 = params["e0"]
    t_arr = params["t_arr"]
    tt_size = params["tt_size"]
    atw_factor = params["atw_factor"]
    f_delta_target = params["f_delta_target"]
    ee_max = params["ee_max"]
    norm = integrate.quad(f_delta_target, e0, ee_max)[0]
    t_size = len(t_arr)
    fc_vec = np.zeros(t_size, dtype=np.float64)
    f_values = {}
    def compute_t(t):
        t = int(t)
        if t not in f_values:
            def f(e):
                return np.exp(- e * t)
            v = integrate.quad(f, e0, ee_max)[0]
            f_values[t] = v
        return f_values[t]
    for t_idx, t in enumerate(t_arr):
        v = compute_t(t)
        if tt_size is not None:
            v = v + atw_factor * compute_t(tt_size - t)
        fc_vec[t_idx] = v
    fc_vec = jnp.array(fc_vec)
    params[tag] = (norm, fc_vec,)
    return norm, fc_vec

def build_hlt_aa_const(params):
    tag = "aa_const"
    if tag in params:
        return params[tag]
    f_e_weight_log = get_f_e_weight_log(params)
    e0 = params["e0"]
    f_delta_target = params["f_delta_target"]
    ee_max = params["ee_max"]
    def f(e):
        return np.exp(f_e_weight_log(e)) * f_delta_target(e)**2
    aa_const = integrate.quad(f, e0, ee_max)[0]
    params[tag] = aa_const
    return aa_const

def aa_from_g(g_t_arr, params):
    """
    g_t_arr.shape == (t_size,)
    delta_target.shape == (n_energies,)
    """
    t_arr = params["t_arr"]
    t_size = len(t_arr)
    aa_mat = build_hlt_aa_mat(params)
    f_vec = build_hlt_f_vec(params)
    aa_const = build_hlt_aa_const(params)
    s = aa_const
    s += (f_vec * g_t_arr).sum()
    s += (aa_mat * g_t_arr * g_t_arr[:, None]).sum()
    return s

def normalization_constraint(g_t_arr, params):
    """
    return new_g_t_arr, constraint_penalty
    return a new g_t_arr which satisfy constraint and penalty reflect the deviation
    """
    does_have_constraint = params["does_have_constraint"]
    if not does_have_constraint:
        return g_t_arr, 0.0
    t_arr = params["t_arr"]
    norm, fc_vec = build_hlt_fc_vec(params)
    g1_t_arr = jnp.ones_like(g_t_arr)
    s_target = norm
    s = (fc_vec * g_t_arr).sum()
    s1 = (fc_vec * g1_t_arr).sum()
    coef = (s_target - s) / s1
    new_g_t_arr = g_t_arr + coef * g1_t_arr
    constraint_penalty = 10.0 * coef**2
    return new_g_t_arr, constraint_penalty

def ww_from_g(g_t_arr, params):
    """
    g_t_arr.shape == (t_size,)
    delta_target.shape == (n_energies,)
    """
    cov = params["cov"]
    t_arr = params["t_arr"]
    new_g_t_arr, constraint_penalty = normalization_constraint(g_t_arr, params)
    g0_t_arr = jnp.zeros(new_g_t_arr.shape, jnp.float64)
    aa = aa_from_g(new_g_t_arr, params) / aa_from_g(g0_t_arr, params)
    ee = get_cov_term(new_g_t_arr, cov) * params["lambda"]
    return aa + ee + constraint_penalty

ww_from_g_wgrad = jax.value_and_grad(ww_from_g)

def mk_g_t_arr_optimization_fcn(params):
    def fcn(g_t_arr, requires_grad=True):
        if requires_grad:
            return ww_from_g_wgrad(g_t_arr, params)
        else:
            return ww_from_g(g_t_arr, params)
    return fcn

@q.timer
def mk_g_t_arr(params):
    fcn = mk_g_t_arr_optimization_fcn(params)
    t_arr = params["t_arr"]
    g_t_arr = params["g_t_arr_init"]
    if g_t_arr is None:
        g_t_arr = np.zeros_like(t_arr, dtype=np.float64)
    for i in range(params["minimization_iter_max"]):
        g_t_arr = q.q_fit_corr.minimize_scipy(fcn, param_arr=g_t_arr)
    return g_t_arr

##############
