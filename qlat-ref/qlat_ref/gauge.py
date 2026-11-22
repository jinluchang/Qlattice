"""
Module ``qlat_ref.gauge``
=========================\n
Pure-NumPy reference implementation of Landau and Coulomb gauge fixing on a
four-dimensional lattice with a plain local iterative scheme (Jacobi
steepest descent with an exponential local gauge transformation; no Fourier
acceleration and no FFT).\n
Documentation: ``docs/qlat-ref/qlat_gauge.md``\n
.. note:: Update the documentation when updating this source file.\n
Array layout and index conventions
----------------------------------\n
A gauge field is a NumPy array ``U`` of shape::\n
    U.shape == (4, Nt, Nx, Ny, Nz, NC, NC)\n
with ``NC = 3``, ``dtype = numpy.complex128`` and (for a link field)\n
    U[mu, t, x, y, z, i, j] = U_mu(t, x, y, z)[i, j]\n
so ``U[mu]`` is the field of links in direction ``mu`` (``0 = t``,
``1 = x``, ``2 = y``, ``3 = z``), its site axes are ``(t, x, y, z)`` and the
last two axes are the (row, column) color indices of the ``SU(3)`` matrix.
A gauge transformation ``g`` is a site field of shape::\n
    g.shape == (Nt, Nx, Ny, Nz, NC, NC)\n
The lattice is periodic in every direction, so a shift by one hop in
direction ``mu`` is ``numpy.roll(..., 1, axis=mu)``: for a link field ``U[mu]``
the coordinate axis of direction ``mu`` is ``axis=mu`` because ``U[mu]``
drops the leading direction axis.\n
The implementation never loops over lattice sites: every site is treated at
once by NumPy array operations on these shapes.  Loops are only over the
four lattice directions and over the (short) Taylor order of the matrix
exponential.
"""

import numpy as np

### -------------------------------------------------------------------

NC = 3

_GELL_MANN = np.asarray(
    [
        [[0, 1, 0], [1, 0, 0], [0, 0, 0]],
        [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]],
        [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
        [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
        [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]],
        [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]],
        [[1, 0, 0], [0, 1, 0], [0, 0, -2]] / np.sqrt(3.0),
    ],
    dtype=np.complex128,
)

### -------------------------------------------------------------------

def get_dirs(dirs=None):
    """
    Return the tuple of link directions to include in the gauge functional.\n
    ``dirs=None`` means all four directions ``(0, 1, 2, 3)`` (Landau gauge);
    ``dirs=(1, 2, 3)`` means the spatial directions only (Coulomb gauge).
    """
    if dirs is None:
        return (0, 1, 2, 3)
    dirs = tuple(int(mu) for mu in dirs)
    assert len(dirs) > 0
    assert len(set(dirs)) == len(dirs)
    for mu in dirs:
        assert 0 <= mu < 4
    return dirs

def get_matrix_adjoint(a):
    """Return the conjugate transpose over the last two (color) axes."""
    return np.swapaxes(np.conjugate(a), -1, -2)

def get_traceless_antihermitian_part(a):
    """
    Return ``(a - a^dagger) / 2`` with its trace removed.\n
    The result is anti-hermitian and traceless; it is the Lie algebra element
    used for the local gauge update.
    """
    h = (a - get_matrix_adjoint(a)) / 2
    tr = np.trace(h, axis1=-2, axis2=-1)
    return h - tr[..., None, None] * np.eye(NC) / NC

def get_su3_exp(h, *, max_order=20):
    """
    Return ``exp(h)`` for an anti-hermitian traceless matrix field ``h``.\n
    The exponential is evaluated with the same truncated Taylor (Horner)
    recurrence as ``qlat_utils.make_matrix_exp`` (truncation order
    ``max_order``), so that the reference and the qlat implementation use the
    identical arithmetic.  For ``exp(h)`` with ``h`` anti-hermitian traceless
    the result is an ``SU(3)`` matrix (up to the Taylor truncation).
    """
    h = np.asarray(h, dtype=np.complex128)
    unit = np.broadcast_to(np.eye(NC), h.shape)
    t2 = h.copy()
    for order in range(max_order, 1, -1):
        t3 = unit + t2 / order
        t2 = h @ t3
    return unit + t2

def get_unit_gauge_field(total_site):
    """Return the unit gauge field ``U_mu(x) = 1`` on ``total_site``."""
    shape = (4,) + tuple(total_site) + (NC, NC)
    U = np.zeros(shape, dtype=np.complex128)
    for mu in range(4):
        U[mu] = np.eye(NC)
    return U

def _get_random_su3_site_field(shape, rng, scale):
    """Return a random ``SU(3)`` site field of the given shape (internal)."""
    coef = rng.standard_normal(tuple(shape) + (8,))
    h = np.zeros(tuple(shape) + (NC, NC), dtype=np.complex128)
    for a in range(8):
        h += 1j * scale * coef[..., a][..., None, None] * (_GELL_MANN[a] / 2)
    return get_su3_exp(h, max_order=30)

def get_random_gauge_field(total_site, *, seed=0, scale=1.0):
    """
    Return a random unitary gauge field on ``total_site``.\n
    The field is deterministic given ``seed``: eight Gaussian random
    coefficients per site and direction are exponentiated into ``SU(3)``
    matrices.  ``scale`` controls the spread of the Gaussian coefficients
    (``scale=1.0`` gives a strongly disordered field, ``scale=0.0`` the unit
    field).
    """
    rng = np.random.default_rng(seed)
    shape = (4,) + tuple(total_site)
    return _get_random_su3_site_field(shape, rng, scale)

def get_random_gauge_transform(total_site, *, seed=0, scale=0.3):
    """
    Return a random ``SU(3)`` gauge transformation on ``total_site``.\n
    The transformation is deterministic given ``seed``.
    """
    rng = np.random.default_rng(seed)
    shape = tuple(total_site)
    return _get_random_su3_site_field(shape, rng, scale)

### -------------------------------------------------------------------

def apply_gauge_transform(U, g):
    """
    Return the gauge-transformed field ``U^g`` with\n
        ``U^g_mu(x) = g(x) U_mu(x) g(x + mu)^dagger``.\n
    ``U`` has shape ``(4, Nt, Nx, Ny, Nz, NC, NC)`` and ``g`` has shape
    ``(Nt, Nx, Ny, Nz, NC, NC)``.
    """
    U = np.asarray(U, dtype=np.complex128)
    g = np.asarray(g, dtype=np.complex128)
    assert U.shape[0] == 4
    assert U.shape[1:] == g.shape
    U_transformed = np.zeros_like(U)
    for mu in range(4):
        g_at_x_plus_mu = np.roll(g, -1, axis=mu)
        U_transformed[mu] = g @ U[mu] @ get_matrix_adjoint(g_at_x_plus_mu)
    return U_transformed

def get_staple_sum(U, *, dirs=None):
    """
    Return the local ``K`` matrix on each site,\n
        ``K(x) = sum_{mu in dirs} (U_mu(x) + U_mu(x - mu)^dagger)``.\n
    It collects the four (or, for Coulomb gauge, three) links that couple to
    the local functional at ``x``.
    """
    U = np.asarray(U, dtype=np.complex128)
    K = np.zeros_like(U[0])
    for mu in get_dirs(dirs):
        K = K + U[mu] + get_matrix_adjoint(np.roll(U[mu], 1, axis=mu))
    return K

def get_local_gauge_gradient(U, *, dirs=None):
    """
    Return the local gradient matrix ``D(x)`` of the gauge functional.\n
    ``D(x)`` is the anti-hermitian traceless part of ``K(x)^dagger`` (with
    ``K`` from :func:`get_staple_sum`); it is the ascent direction of the
    functional and ``exp(alpha * D(x))`` with ``alpha > 0`` is the local gauge
    transformation that increases the functional.  The
    gauge condition is ``D(x) = 0`` for all ``x``.
    """
    K = get_staple_sum(U, dirs=dirs)
    return get_traceless_antihermitian_part(get_matrix_adjoint(K))

def get_gauge_functional(U, *, dirs=None):
    """
    Return the (normalized) gauge functional,\n
        ``F = sum_{x, mu in dirs} Re Tr U_mu(x) / (V * len(dirs) * NC)``,\n
    where ``V`` is the total number of sites.  Landau gauge uses
    ``dirs=(0, 1, 2, 3)``; the Coulomb functional uses ``dirs=(1, 2, 3)``.
    """
    U = np.asarray(U, dtype=np.complex128)
    dirs = get_dirs(dirs)
    volume = int(np.prod(U.shape[1:5]))
    total = 0.0
    for mu in dirs:
        total += np.sum(np.trace(U[mu], axis1=-2, axis2=-1).real)
    return float(total / (volume * len(dirs) * NC))

def get_coulomb_functional_list(U):
    """
    Return the Coulomb functional of each time slice.\n
    The Coulomb functional is a sum over the spatial links only, so it splits
    into one independent ``3``-dimensional Landau functional per time slice.
    """
    U = np.asarray(U, dtype=np.complex128)
    n_t = U.shape[1]
    n_spatial = int(np.prod(U.shape[2:5]))
    total = np.zeros(n_t, dtype=np.float64)
    for mu in (1, 2, 3):
        total += np.sum(np.trace(U[mu], axis1=-2, axis2=-1).real, axis=(1, 2, 3))
    return [float(value) for value in total / (n_spatial * 3 * NC)]

def get_divergence_from_gradient(D):
    """Return the RMS divergence ``theta`` for a local gradient field ``D``."""
    D = np.asarray(D, dtype=np.complex128)
    volume = int(np.prod(D.shape[:4]))
    quad = np.trace(D @ get_matrix_adjoint(D), axis1=-2, axis2=-1).real
    return float(np.sqrt(np.sum(quad) / (volume * NC)))

def get_gauge_divergence(U, *, dirs=None):
    """
    Return the RMS lattice divergence ``theta`` of the gauge field.\n
    ``theta`` is ``sqrt(sum_x Tr[D(x) D(x)^dagger] / (V * NC))`` with ``D`` the
    local gradient of the functional; it vanishes exactly when the gauge
    condition ``D(x) = 0`` holds.  This is the convergence measure of
    :func:`fix_landau_gauge` and :func:`fix_coulomb_gauge`.
    """
    D = get_local_gauge_gradient(U, dirs=dirs)
    return get_divergence_from_gradient(D)

def get_coulomb_divergence_list(U):
    """
    Return the RMS divergence of each time slice for the Coulomb gauge.\n
    Each time slice is an independent ``3``-dimensional Landau problem, so the
    per-slice values are the natural convergence monitor of Coulomb gauge
    fixing.
    """
    D = get_local_gauge_gradient(U, dirs=(1, 2, 3))
    D.shape[0]
    n_spatial = int(np.prod(D.shape[1:4]))
    quad = np.trace(D @ get_matrix_adjoint(D), axis1=-2, axis2=-1).real
    per_slice = np.sum(quad, axis=(1, 2, 3))
    return [float(value) for value in np.sqrt(per_slice / (n_spatial * NC))]

def get_unitarity_error(U):
    """Return the maximal ``|U U^dagger - 1|`` matrix element of the field."""
    U = np.asarray(U, dtype=np.complex128)
    deviation = U @ get_matrix_adjoint(U) - np.eye(NC)
    return float(np.max(np.abs(deviation)))

def check_unitarity(U, *, eps=1e-12):
    """Return whether every link of ``U`` is unitary within ``eps``."""
    return get_unitarity_error(U) <= eps

### -------------------------------------------------------------------

def _fix_gauge(U, *, dirs, gauge_name, max_iter, eps, alpha, verbose):
    """
    Shared iteration core of Landau and Coulomb gauge fixing (internal).\n
    The local gauge transformation of one sweep is\n
        ``g(x) = exp(alpha * D(x))``   (``D`` from
        :func:`get_local_gauge_gradient`),\n
    applied to every site at once (Jacobi sweep):\n
        ``U_mu(x) <- g(x) U_mu(x) g(x + mu)^dagger``.\n
    The iteration stops when the RMS divergence is at most ``eps`` or after
    ``max_iter`` sweeps.
    """
    U = np.array(U, dtype=np.complex128, copy=True)
    assert alpha > 0.0
    dirs = get_dirs(dirs)
    functional_list = []
    divergence_list = []
    converged = False
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        D = get_local_gauge_gradient(U, dirs=dirs)
        divergence = get_divergence_from_gradient(D)
        functional = get_gauge_functional(U, dirs=dirs)
        functional_list.append(functional)
        divergence_list.append(divergence)
        if verbose:
            print(
                f"qlat_ref: {gauge_name} sweep {n_iter}:"
                f" functional = {functional:.15e}"
                f" divergence = {divergence:.3e}"
            )
        if divergence <= eps:
            converged = True
            break
        g = get_su3_exp(alpha * D)
        U = apply_gauge_transform(U, g)
    D = get_local_gauge_gradient(U, dirs=dirs)
    info = {
        "gauge": gauge_name,
        "dirs": dirs,
        "converged": converged,
        "n_iter": n_iter,
        "max_iter": max_iter,
        "eps": eps,
        "alpha": alpha,
        "functional": get_gauge_functional(U, dirs=dirs),
        "divergence": get_divergence_from_gradient(D),
        "functional_list": functional_list,
        "divergence_list": divergence_list,
    }
    return U, info

def fix_landau_gauge(U, *, max_iter=10000, eps=1e-10, alpha=0.1, verbose=False):
    """
    Return ``(U_fixed, info)`` with ``U_fixed`` in Landau gauge.\n
    Landau gauge maximizes the functional\n
        ``F = (1 / (V * 4 * NC)) sum_{x, mu} Re Tr U_mu(x)``\n
    over all gauge transformations; the plain local iteration of
    :func:`_fix_gauge` uses all four directions.  The iteration stops when the
    RMS divergence is at most ``eps`` (or after ``max_iter`` sweeps); the
    returned ``info`` dictionary records the convergence history.
    """
    return _fix_gauge(
        U,
        dirs=(0, 1, 2, 3),
        gauge_name="landau",
        max_iter=max_iter,
        eps=eps,
        alpha=alpha,
        verbose=verbose,
    )

def fix_coulomb_gauge(U, *, max_iter=10000, eps=1e-10, alpha=0.1, verbose=False):
    """
    Return ``(U_fixed, info)`` with ``U_fixed`` in Coulomb gauge.\n
    Coulomb gauge maximizes the functional\n
        ``F_C = (1 / (V * 3 * NC)) sum_{x, i} Re Tr U_i(x)``\n
    over gauge transformations ``g(t, x)`` that act independently on every
    time slice; the iteration therefore uses the spatial directions only
    (``i = 1, 2, 3``) and is a separate 3D Landau fixing per time slice.  The
    gauge transformation is still applied to all four link directions, so the
    result is a complete 4D gauge transformation.  The iteration stops when
    the RMS divergence is at most ``eps`` (or after ``max_iter`` sweeps).
    """
    return _fix_gauge(
        U,
        dirs=(1, 2, 3),
        gauge_name="coulomb",
        max_iter=max_iter,
        eps=eps,
        alpha=alpha,
        verbose=verbose,
    )
