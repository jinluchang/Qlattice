# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.gauge_fixing``
============================\n
Plain local steepest-descent **Landau** and **Coulomb** gauge fixing on a
distributed ``GaugeField`` (``Field<ColorMatrix>``): no Fourier acceleration
and no FFT.\n
The gauge functional is\n
    F = sum_{x, mu in dirs} Re Tr U_mu(x) / (V * len(dirs) * 3)\n
and the local update of one Jacobi sweep is ``g(x) = exp(alpha * D(x))`` with
``D(x)`` the anti-hermitian traceless part of the staple-sum matrix
``K(x) = sum_{mu in dirs} (U_mu(x) + U_mu(x-mu)^dagger)``, followed by\n
    U_mu(x) -> g(x) U_mu(x) g(x+mu)^dagger.\n
The iteration stops when the RMS divergence ``theta`` is at most ``eps`` or
after ``max_iter`` sweeps.\n
Documentation: ``docs/qlat/qlat_gauge_fixing.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from .qcd cimport GaugeField

### -------------------------------------------------------------------

cdef dict _get_info_dict(GaugeFixingInfo info, str gauge_name):
    return {
        "gauge": gauge_name,
        "converged": info.converged == True,
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
    Int max_iter=10000,
    RealD eps=1.0e-10,
    RealD alpha=0.1,
    Bool verbose=False,
):
    """
    Plain local steepest-descent gauge fixing of ``gf`` in place with the link
    directions ``dirs``.\n
    Returns a dictionary with the convergence history: ``gauge``,
    ``converged``, ``n_iter``, ``max_iter``, ``eps``, ``alpha``,
    ``functional``, ``divergence``, ``functional_list`` and
    ``divergence_list``.
    """
    cdef std_vector[Int] c_dirs
    for mu in dirs:
        c_dirs.push_back(mu)
    cdef GaugeFixingInfo info = gf_fix_gauge(
        gf.xxx().val(), c_dirs, max_iter, eps, alpha, verbose
    )
    return _get_info_dict(info, "custom")

def fix_landau_gauge(
    GaugeField gf,
    Int max_iter=10000,
    RealD eps=1.0e-10,
    RealD alpha=0.1,
    Bool verbose=False,
):
    """
    Plain local steepest-descent Landau gauge fixing of ``gf`` in place.\n
    Maximizes the normalized functional summed over all four directions.  The
    iteration stops when the RMS divergence is at most ``eps`` or after
    ``max_iter`` sweeps.  Returns a dictionary with the convergence history.
    """
    cdef GaugeFixingInfo info = gf_fix_landau_gauge(
        gf.xxx().val(), max_iter, eps, alpha, verbose
    )
    return _get_info_dict(info, "landau")

def fix_coulomb_gauge(
    GaugeField gf,
    Int max_iter=10000,
    RealD eps=1.0e-10,
    RealD alpha=0.1,
    Bool verbose=False,
):
    """
    Plain local steepest-descent Coulomb gauge fixing of ``gf`` in place.\n
    Maximizes the normalized functional summed over the spatial directions,
    which splits into independent 3D Landau fixings per time slice; the
    resulting gauge transformation is still applied to all four link
    directions.  The iteration stops when the RMS divergence is at most
    ``eps`` or after ``max_iter`` sweeps.  Returns a dictionary with the
    convergence history.
    """
    cdef GaugeFixingInfo info = gf_fix_coulomb_gauge(
        gf.xxx().val(), max_iter, eps, alpha, verbose
    )
    return _get_info_dict(info, "coulomb")

def get_gauge_functional(GaugeField gf, dirs):
    """
    Return the normalized gauge functional
    ``sum_{x, mu in dirs} Re Tr U_mu(x) / (V * len(dirs) * 3)``.
    """
    cdef std_vector[Int] c_dirs
    for mu in dirs:
        c_dirs.push_back(mu)
    return gf_gauge_functional(gf.xxx().val(), c_dirs)

def get_gauge_divergence(GaugeField gf, dirs):
    """
    Return the RMS lattice divergence ``theta`` of the gauge field for the
    link directions ``dirs``; it vanishes exactly in the fixed gauge.
    """
    cdef std_vector[Int] c_dirs
    for mu in dirs:
        c_dirs.push_back(mu)
    return gf_gauge_divergence(gf.xxx().val(), c_dirs)

def get_unitarity_error(GaugeField gf):
    """Return the maximal ``|U U^dagger - 1|`` matrix element of the field."""
    return gf_gauge_unitarity_error(gf.xxx().val())
