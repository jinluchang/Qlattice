__all__ = [
        'mk_hlt_params',
        #
        'delta_from_g',
        #
        'aa_from_g_via_sum',
        'ww_from_g_via_sum',
        'ww_from_g_wgrad_via_sum',
        'mk_g_t_arr_via_sum',
        #
        'aa_from_g',
        'ww_from_g',
        'ww_from_g_wgrad',
        'mk_g_t_arr',
        #
        'get_f_e_weight',
        ]

import jax
import jax.numpy as jnp

import numpy as np
import scipy.integrate as integrate

import qlat_utils as q

if jnp.zeros(2, dtype=jnp.float64).dtype != jnp.float64:
    raise Exception(f"ERROR: double precision not available in JAX. Please set: 'export JAX_ENABLE_X64=True' to use double precision number here.")

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
    params["lambda"] = 1.0
    params["alpha"] = -0.01 # override by "f_e_weight", see `get_f_e_weight`
    params["f_e_weight"] = None
    params["tt_size"] = None
    params["atw_factor"] = 1.0 # only has effects if "tt_size" is not None
    params["mini_iter_max"] = 1
    params["g_t_arr_init"] = None
    return params

def get_f_e_weight(params):
    f = params["f_e_weight"]
    if f is not None:
        return f
    alpha = params["alpha"]
    def f(e):
        return np.exp(alpha * e)
    return f

@jax.jit
def delta_from_g(g_t_arr, t_arr, e_arr):
    """
    return delta
    g_t_arr.shape == (..., t_size,)
    delta.shape[-1] == (..., n_energies,)
    """
    delta = (jnp.exp(-e_arr[:, None] * t_arr) * g_t_arr[..., None, :]).sum(-1)
    return delta

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
    f_e_weight = get_f_e_weight(params)
    e0 = params["e0"]
    e_w = f_e_weight(e_arr)
    e_w = jnp.where(e_arr >= e0, e_w, 0.0)
    delta = delta_from_g(g_t_arr, t_arr, e_arr)
    if tt_size is not None:
        delta = delta + atw_factor * delta_from_g(g_t_arr, tt_size - t_arr, e_arr)
    f_e = (delta - delta_target)**2 * e_w
    s1 = f_e[0] * (e_arr[1] - e_arr[0]) / 2
    s2 = f_e[-1] * (e_arr[-1] - e_arr[-2]) / 2
    ss = f_e[1:-1] * (e_arr[2:] - e_arr[:-2]) / 2
    return s1 + s2 + ss.sum()

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

def ww_from_g_via_sum(g_t_arr, params):
    """
    g_t_arr.shape == (t_size,)
    delta_target.shape == (n_energies,)
    """
    cov = params["cov"]
    t_arr = params["t_arr"]
    t_size = len(t_arr)
    aa = aa_from_g_via_sum(g_t_arr, params)
    ee = get_cov_term(g_t_arr, cov) * params["lambda"]
    return aa + ee

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
    for i in range(params["mini_iter_max"]):
        g_t_arr = q.q_fit_corr.minimize_scipy(fcn, param_arr=g_t_arr)
    return g_t_arr

############## with integration ##############

def build_hlt_aa_mat(params):
    tag = "aa_mat"
    if tag in params:
        return params[tag]
    f_e_weight = get_f_e_weight(params)
    e0 = params["e0"]
    t_arr = params["t_arr"]
    tt_size = params["tt_size"]
    atw_factor = params["atw_factor"]
    t_size = len(t_arr)
    aa_mat = np.zeros((t_size, t_size,), dtype=np.float64)
    aa_values = {}
    def compute_t_sum(t_sum):
        t_sum = int(t_sum)
        if t_sum not in aa_values:
            def f(e):
                return f_e_weight(e) * np.exp(-e * t_sum)
            v = integrate.quad(f, e0, np.inf)[0]
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
    f_e_weight = get_f_e_weight(params)
    e0 = params["e0"]
    t_arr = params["t_arr"]
    tt_size = params["tt_size"]
    atw_factor = params["atw_factor"]
    f_delta_target = params["f_delta_target"]
    t_size = len(t_arr)
    f_vec = np.zeros(t_size, dtype=np.float64)
    f_values = {}
    def compute_t(t):
        t = int(t)
        if t not in f_values:
            def f(e):
                return -2 * f_e_weight(e) * f_delta_target(e) * np.exp(-e * t)
            v = integrate.quad(f, e0, np.inf)[0]
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

def build_hlt_aa_const(params):
    tag = "aa_const"
    if tag in params:
        return params[tag]
    f_e_weight = get_f_e_weight(params)
    e0 = params["e0"]
    f_delta_target = params["f_delta_target"]
    def f(e):
        return f_e_weight(e) * f_delta_target(e)**2
    aa_const = integrate.quad(f, e0, np.inf)[0]
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

def ww_from_g(g_t_arr, params):
    """
    g_t_arr.shape == (t_size,)
    delta_target.shape == (n_energies,)
    """
    cov = params["cov"]
    t_arr = params["t_arr"]
    t_size = len(t_arr)
    aa = aa_from_g(g_t_arr, params)
    ee = get_cov_term(g_t_arr, cov) * params["lambda"]
    return aa + ee

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
    for i in range(params["mini_iter_max"]):
        g_t_arr = q.q_fit_corr.minimize_scipy(fcn, param_arr=g_t_arr)
    return g_t_arr
