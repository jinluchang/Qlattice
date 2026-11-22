# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat.all cimport *
from . cimport everything as cc
from qlat.qcd cimport GaugeField

### -------------------------------------------------------------------

cdef dict _get_info_dict(cc.GaugeFixingInfo info, str gauge_name):
    return {
        "gauge": gauge_name,
        "converged": bool(info.converged),
        "n_iter": int(info.n_iter),
        "max_iter": int(info.max_iter),
        "eps": float(info.eps),
        "alpha": float(info.alpha),
        "functional": float(info.functional),
        "divergence": float(info.divergence),
        "functional_list": [float(x) for x in info.functional_list],
        "divergence_list": [float(x) for x in info.divergence_list],
    }

### -------------------------------------------------------------------

def fix_gauge(
    GaugeField gf,
    dirs,
    cc.Int max_iter=10000,
    cc.RealD eps=1.0e-10,
    cc.RealD alpha=0.1,
    cc.Bool verbose=False,
):
    """
    Plain local steepest-descent gauge fixing of ``gf`` in place with the
    link directions ``dirs``.\n
    The local update is ``g(x) = exp(alpha * D(x))`` with ``D`` the
    anti-hermitian traceless part of the staple-sum matrix; the iteration
    stops when the RMS divergence is at most ``eps`` or after ``max_iter``
    sweeps.  Returns a dictionary with the convergence history.
    """
    cdef cc.std_vector[cc.Int] c_dirs
    for mu in dirs:
        c_dirs.push_back(mu)
    cdef cc.GaugeFixingInfo info = cc.gf_fix_gauge(
        gf.xxx().val(), c_dirs, max_iter, eps, alpha, verbose
    )
    return _get_info_dict(info, "custom")

def fix_landau_gauge(
    GaugeField gf,
    cc.Int max_iter=10000,
    cc.RealD eps=1.0e-10,
    cc.RealD alpha=0.1,
    cc.Bool verbose=False,
):
    """
    Plain local steepest-descent Landau gauge fixing of ``gf`` in place.\n
    Maximizes ``sum_{x, mu} Re Tr U_mu(x)`` over all four directions.  The
    iteration stops when the RMS divergence is at most ``eps`` or after
    ``max_iter`` sweeps.  Returns a dictionary with the convergence history.
    """
    cdef cc.GaugeFixingInfo info = cc.gf_fix_landau_gauge(
        gf.xxx().val(), max_iter, eps, alpha, verbose
    )
    return _get_info_dict(info, "landau")

def fix_coulomb_gauge(
    GaugeField gf,
    cc.Int max_iter=10000,
    cc.RealD eps=1.0e-10,
    cc.RealD alpha=0.1,
    cc.Bool verbose=False,
):
    """
    Plain local steepest-descent Coulomb gauge fixing of ``gf`` in place.\n
    Maximizes ``sum_{x, i} Re Tr U_i(x)`` over the spatial directions, which
    splits into independent 3D Landau fixings per time slice.  The iteration
    stops when the RMS divergence is at most ``eps`` or after ``max_iter``
    sweeps.  Returns a dictionary with the convergence history.
    """
    cdef cc.GaugeFixingInfo info = cc.gf_fix_coulomb_gauge(
        gf.xxx().val(), max_iter, eps, alpha, verbose
    )
    return _get_info_dict(info, "coulomb")

### -------------------------------------------------------------------

def get_gauge_functional(GaugeField gf, dirs):
    """
    Return the normalized gauge functional
    ``sum_{x, mu in dirs} Re Tr U_mu(x) / (V * len(dirs) * 3)``.
    """
    cdef cc.std_vector[cc.Int] c_dirs
    for mu in dirs:
        c_dirs.push_back(mu)
    return cc.gf_gauge_functional(gf.xxx().val(), c_dirs)

def get_gauge_divergence(GaugeField gf, dirs):
    """
    Return the RMS lattice divergence ``theta`` of the gauge field for the
    link directions ``dirs``.
    """
    cdef cc.std_vector[cc.Int] c_dirs
    for mu in dirs:
        c_dirs.push_back(mu)
    return cc.gf_gauge_divergence(gf.xxx().val(), c_dirs)

def get_unitarity_error(GaugeField gf):
    """Return the maximal ``|U U^dagger - 1|`` matrix element of the field."""
    return cc.gf_gauge_unitarity_error(gf.xxx().val())
