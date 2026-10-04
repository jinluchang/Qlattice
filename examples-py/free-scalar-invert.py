#!/usr/bin/env python3

# Tests for q.free_scalar_invert, q.free_scalar_mom_invert,
# q.free_scalar_deriv_mom and q.free_scalar_invert_deriv.
#
# free_scalar_deriv_mom multiplies a momentum-space field by
#     prod_mu [ d_even(k_mu)^m (i sin(k_mu))^e ],  deriv[mu] = 2 m + e,
# i.e. the leftover odd order is one symmetric (central) difference factor
# i sin(k_mu), while the even kernel d_even is selected by even_deriv_kernel:
# "half" (the default) is the laplacian power -4 sin^2(k_mu/2) and "central"
# is the square of the odd factor, -sin^2(k_mu).
#
# free_scalar_invert(src, mass) must solve the free lattice scalar
# equation
#     (4 sinh^2(mass/2) - laplacian) sol = src
# on the periodic lattice, i.e. multiply src~(k) by
#     1 / (4 sinh^2(mass/2) + sum_mu 4 sin^2(k_mu/2))
# in momentum space.  The checks below pin down
#   * the pole-mass convention 4 sinh^2(mass/2) (and not mass^2),
#   * the symmetric 1/sqrt(V) normalization of the forward/inverse FFT pair,
#   * the use of global (not local) lattice coordinates as the momentum label,
#     which is what makes the function correct in multi-node runs.
# An independent numpy DFT reference is used, so a regression in the qlat FFT
# implementation is caught here as well.

import numpy as np

import qlat as q

check_eps = 1e-10

total_site = q.Coordinate([4, 4, 4, 8])
mass = 0.3

# The test runner uses `mpiexec -n 2` and begin_with_mpi() takes the first entry
# with the matching node count: split spatially so that the global coordinate
# of a site differs from its local coordinate.
size_node_list = [
    [2, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 1],
    [2, 2, 2, 2],
]

def get_latt_size(geo):
    return tuple(int(geo.total_site[i]) for i in range(4))

def get_m_pi_sq(mass):
    # D(k) = 4 sinh^2(mass/2) + sum_mu 4 sin^2(k_mu/2) = 2 (cosh(mass) - cos k)
    return 4.0 * np.sinh(mass / 2.0) ** 2

def mk_s2_grid(latt_size):
    # sum_mu 4 sin^2(k_mu/2) on the momentum-index grid n_mu = 0 .. L_mu - 1
    ax = []
    for i in range(4):
        k = 2.0 * np.pi * np.arange(latt_size[i]) / latt_size[i]
        ax.append(4.0 * np.sin(k / 2.0) ** 2)
    return (
        ax[0][:, None, None, None]
        + ax[1][None, :, None, None]
        + ax[2][None, None, :, None]
        + ax[3][None, None, None, :]
    )

def mk_mom_factor(latt_size, mass):
    return get_m_pi_sq(mass) + mk_s2_grid(latt_size)

def mk_mom_factor_twist(latt_size, mass, momtwist):
    # same factor with the momentum grid shifted by the twist angles
    ax = []
    for i in range(4):
        k = 2.0 * np.pi * (np.arange(latt_size[i]) + momtwist[i]) / latt_size[i]
        ax.append(4.0 * np.sin(k / 2.0) ** 2)
    s2 = (
        ax[0][:, None, None, None]
        + ax[1][None, :, None, None]
        + ax[2][None, None, :, None]
        + ax[3][None, None, None, :]
    )
    return get_m_pi_sq(mass) + s2

def ref_free_scalar_invert(src_arr, latt_size, mass):
    # independent DFT reference: (1/V) sum_k e^{i k x} src~(k) / D(k)
    return np.fft.ifftn(np.fft.fftn(src_arr) / mk_mom_factor(latt_size, mass))

def mk_global_src(rs, latt_size):
    # deterministic source defined by global coordinates, so that the same
    # global field is used for any node layout
    re = rs.split("re").u_rand_arr(latt_size) - 0.5
    im = rs.split("im").u_rand_arr(latt_size) - 0.5
    return re + 1j * im

def mk_delta_src(latt_size):
    arr = np.zeros(latt_size, dtype=np.complex128)
    arr[0, 0, 0, 0] = 1.0
    return arr

def mk_field_from_global(geo, arr):
    f = q.FieldComplexD(geo, 1)
    buf = np.asarray(f)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        buf[index, 0] = arr[xg[0], xg[1], xg[2], xg[3]]
    return f

def get_global_arr(geo, f):
    # gather local data into a global array indexed by global coordinates;
    # applied to a momentum-space field it is indexed by momentum
    buf = np.asarray(f)
    arr = np.zeros(get_latt_size(geo), dtype=np.complex128)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        arr[xg[0], xg[1], xg[2], xg[3]] = buf[index, 0]
    # each global site is owned by exactly one node, so the sum is the gather
    return q.glb_sum(arr)

def get_neg_lap(arr):
    # -laplacian = sum_mu (2 - shift_mu - shift_mu^-1)
    s = np.zeros_like(arr)
    for i in range(4):
        s += 2.0 * arr - np.roll(arr, 1, axis=i) - np.roll(arr, -1, axis=i)
    return s

def mk_signed_index(latt_size, i):
    # smod(n_i, L_i): the folded momentum index used by the C++ kernels
    L = latt_size[i]
    n = np.arange(L)
    return ((n + L // 2) % L) - L // 2

def mk_k_axis(latt_size, momtwist, i):
    # k_i = 2 pi ( smod(n_i, L_i) + momtwist_i ) / L_i
    return 2.0 * np.pi * (mk_signed_index(latt_size, i) + momtwist[i]) / latt_size[i]

def mk_deriv_factor(latt_size, deriv, momtwist, even_deriv_kernel=None):
    # independent reference: prod_mu [ d_even(k_mu)^m (i sin(k_mu))^e ]
    # with deriv[mu] = 2 m + e and the odd (central difference) factor
    # i sin(k_mu), which vanishes at the self-conjugate momentum k_mu = pi.
    # even_deriv_kernel selects d_even: "half" (default) is the laplacian power
    # -4 sin^2(k_mu/2), "central" is the square of the odd factor -sin^2(k_mu)
    if even_deriv_kernel is None:
        even_deriv_kernel = "half"
    fac = np.ones(latt_size, dtype=np.complex128)
    for i in range(4):
        n_signed = mk_signed_index(latt_size, i)
        k = mk_k_axis(latt_size, momtwist, i)
        shape = [1, 1, 1, 1]
        shape[i] = latt_size[i]
        is_edge = np.abs(2.0 * (n_signed + momtwist[i])) == latt_size[i]
        if even_deriv_kernel == "half":
            d_even = (-4.0 * np.sin(k / 2.0) ** 2).reshape(shape)
        elif even_deriv_kernel == "central":
            d_even = np.where(is_edge, 0.0, -(np.sin(k) ** 2)).reshape(shape)
        else:
            raise Exception(
                f"mk_deriv_factor: bad even_deriv_kernel={even_deriv_kernel}"
            )
        d_odd = np.where(is_edge, 0.0, 1.0j * np.sin(k)).reshape(shape)
        n = int(deriv[i])
        for _ in range(n // 2):
            fac = fac * d_even
        if n % 2 == 1:
            fac = fac * d_odd
    return fac

def mk_self_conjugate_mask(latt_size, mu, momtwist):
    # True at the self-conjugate momentum k_mu = pi, where the odd (central
    # difference) factor i sin(k_mu) vanishes
    edge = (
        np.abs(2.0 * (mk_signed_index(latt_size, mu) + momtwist[mu])) == latt_size[mu]
    )
    shape = [1, 1, 1, 1]
    shape[mu] = latt_size[mu]
    return np.broadcast_to(edge.reshape(shape), latt_size)

q.begin_with_mpi(size_node_list)

geo = q.Geometry(total_site)
latt_size = get_latt_size(geo)
q.json_results_append(f"free-scalar-invert: geo.show()={geo.show()}")
rs = q.RngState("free-scalar-invert")

# --- delta source: compare with an independent numpy DFT reference
delta_src = mk_field_from_global(geo, mk_delta_src(latt_size))
delta_sol = q.free_scalar_invert(delta_src, mass)
delta_sol_arr = get_global_arr(geo, delta_sol)
delta_ref_arr = ref_free_scalar_invert(mk_delta_src(latt_size), latt_size, mass)
err_delta = float(np.max(np.abs(delta_sol_arr - delta_ref_arr)))
assert err_delta < check_eps, err_delta
q.json_results_append(
    f"free-scalar-invert: delta matches DFT reference = {err_delta < check_eps}"
)
q.json_results_append(
    "free-scalar-invert: delta sol(xg=0)",
    float(delta_sol_arr[0, 0, 0, 0].real),
    check_eps,
)

# --- pole mass: the spatial zero mode of a delta source is proportional to
# cosh(mass (t - total_site[3] / 2)), so G(3)/G(4) == cosh(mass) exactly
g_t = delta_sol_arr.sum(axis=(0, 1, 2))
pole_ratio = float(np.real(g_t[3] / g_t[4]))
assert abs(pole_ratio - np.cosh(mass)) < check_eps, pole_ratio
q.json_results_append("free-scalar-invert: delta G(3)/G(4)", pole_ratio, check_eps)

# --- random source: DFT reference, operator residual, src is not modified
src_arr = mk_global_src(rs.split("src"), latt_size)
src = mk_field_from_global(geo, src_arr)
src_before = get_global_arr(geo, src)
sol = q.free_scalar_invert(src, mass)
sol_arr = get_global_arr(geo, sol)
ref_arr = ref_free_scalar_invert(src_arr, latt_size, mass)
err_ref = float(np.max(np.abs(sol_arr - ref_arr)))
assert err_ref < check_eps, err_ref
err_src = float(np.max(np.abs(get_global_arr(geo, src) - src_before)))
assert err_src == 0.0, err_src
q.json_results_append(
    f"free-scalar-invert: rand matches DFT reference = {err_ref < check_eps}"
)
q.json_results_append(
    f"free-scalar-invert: rand leaves src unchanged = {err_src == 0.0}"
)

# --- (4 sinh^2(mass/2) - laplacian) sol = src, i.e. the mass convention
res_arr = get_m_pi_sq(mass) * sol_arr + get_neg_lap(sol_arr) - src_arr
err_op = float(np.max(np.abs(res_arr)))
assert err_op < check_eps, err_op
q.json_results_append(f"free-scalar-invert: operator residual = {err_op < check_eps}")
q.json_results_append(
    "free-scalar-invert: rand sol(xg=(1,2,3,4))",
    np.array([sol_arr[1, 2, 3, 4].real, sol_arr[1, 2, 3, 4].imag]),
    check_eps,
)

# --- mode_fft=0 must agree with the default mode_fft=1
sol0 = q.free_scalar_invert(src, mass, mode_fft=0)
err_mode = float(np.max(np.abs(get_global_arr(geo, sol0) - sol_arr)))
assert err_mode < check_eps, err_mode
q.json_results_append(
    f"free-scalar-invert: mode_fft=0 vs mode_fft=1 = {err_mode < check_eps}"
)

# --- free_scalar_mom_invert applies 1 / D(k) in momentum space
fft_f = q.mk_fft(is_forward=True, is_normalizing=True)
fft_b = q.mk_fft(is_forward=False, is_normalizing=True)
mom = fft_f * src
q.free_scalar_mom_invert(mom, mass)
mom_arr = get_global_arr(geo, mom)
mom_ref_arr = (
    np.fft.fftn(src_arr) / np.sqrt(np.prod(latt_size)) / mk_mom_factor(latt_size, mass)
)
err_mom = float(np.max(np.abs(mom_arr - mom_ref_arr)))
assert err_mom < check_eps, err_mom
q.json_results_append(
    f"free-scalar-invert: mom matches DFT reference = {err_mom < check_eps}"
)
q.json_results_append(
    "free-scalar-invert: mom sol(k=(1,2,3,4))",
    np.array([mom_arr[1, 2, 3, 4].real, mom_arr[1, 2, 3, 4].imag]),
    check_eps,
)

# --- momtwist is exposed and shifts the momentum grid used by the kernel
momtwist_list = [0.5, -0.25, 0.0, 0.125]
momtwist = q.CoordinateD(momtwist_list)
mom_tw = fft_f * src
q.free_scalar_mom_invert(mom_tw, mass, momtwist)
mom_tw_arr = get_global_arr(geo, mom_tw)
mom_tw_ref_arr = (
    np.fft.fftn(src_arr)
    / np.sqrt(np.prod(latt_size))
    / mk_mom_factor_twist(latt_size, mass, momtwist_list)
)
err_mom_tw = float(np.max(np.abs(mom_tw_arr - mom_tw_ref_arr)))
assert err_mom_tw < check_eps, err_mom_tw
q.json_results_append(
    f"free-scalar-invert: mom momtwist matches DFT reference = {err_mom_tw < check_eps}"
)
q.json_results_append(
    "free-scalar-invert: mom momtwist sol(k=(1,2,3,4))",
    np.array([mom_tw_arr[1, 2, 3, 4].real, mom_tw_arr[1, 2, 3, 4].imag]),
    check_eps,
)

# --- free_scalar_deriv_mom: the bare lattice derivative factor in momentum
#     space.  Each order deriv[mu] = 2 m + e is split into an even power of
#     d_even(k_mu) and, for the leftover odd order, one symmetric (central)
#     difference factor i sin(k_mu).  even_deriv_kernel selects d_even:
#     "half" (the default) is the laplacian power -4 sin^2(k_mu/2), so
#     deriv[mu] = 2 is minus the mu term of the lattice laplacian, while
#     "central" is the square of the odd factor, -sin^2(k_mu).  It is checked
#     against (i) an independent momentum-space DFT reference, (ii) the local
#     central difference (and its square) in position space, and (iii) the
#     second derivative.
zero_twist = [0.0, 0.0, 0.0, 0.0]
src_mom_arr = np.fft.fftn(src_arr) / np.sqrt(np.prod(latt_size))

# deriv=None (all orders zero) must leave the field untouched
f_id = fft_f * src
id_before = get_global_arr(geo, f_id)
q.free_scalar_deriv_mom(f_id, None)
err_id = float(np.max(np.abs(get_global_arr(geo, f_id) - id_before)))
assert err_id == 0.0, err_id
q.json_results_append(
    f"free-scalar-invert: deriv none leaves field unchanged = {err_id == 0.0}"
)

# momentum-space DFT reference, single and mixed directions, with and without
# a nonzero momentum twist
deriv_cases = [
    ("dx", [1, 0, 0, 0], zero_twist),
    ("dt", [0, 0, 0, 1], zero_twist),
    ("dx dy", [1, 1, 0, 0], zero_twist),
    ("dx twist", [1, 0, 0, 0], momtwist_list),
    ("dx3", [3, 0, 0, 0], zero_twist),
]
deriv_out = {}
for tag, deriv, twist_list in deriv_cases:
    f_mom = fft_f * src
    q.free_scalar_deriv_mom(f_mom, deriv, q.CoordinateD(twist_list))
    out_arr = get_global_arr(geo, f_mom)
    deriv_out[tag] = out_arr
    ref_deriv_arr = src_mom_arr * mk_deriv_factor(latt_size, deriv, twist_list)
    err_deriv = float(np.max(np.abs(out_arr - ref_deriv_arr)))
    assert err_deriv < check_eps, (tag, err_deriv)
    q.json_results_append(
        f"free-scalar-invert: deriv {tag} matches DFT reference = {err_deriv < check_eps}"
    )
q.json_results_append(
    "free-scalar-invert: deriv dx sol(k=(1,2,3,4))",
    np.array([deriv_out["dx"][1, 2, 3, 4].real, deriv_out["dx"][1, 2, 3, 4].imag]),
    check_eps,
)

# even_deriv_kernel="central" replaces the even laplacian power by the square of
# the odd central difference factor; check it against the same DFT reference
central_cases = [
    ("d0^2 central", [2, 0, 0, 0]),
    ("d3^2 central", [0, 0, 0, 2]),
    ("dx3 central", [3, 0, 0, 0]),
    ("d2 dx central", [2, 1, 0, 0]),
]
for tag, deriv in central_cases:
    f_mom = fft_f * src
    q.free_scalar_deriv_mom(f_mom, deriv, even_deriv_kernel="central")
    out_arr = get_global_arr(geo, f_mom)
    ref_deriv_arr = src_mom_arr * mk_deriv_factor(
        latt_size, deriv, zero_twist, "central"
    )
    err_deriv = float(np.max(np.abs(out_arr - ref_deriv_arr)))
    assert err_deriv < check_eps, (tag, err_deriv)
    q.json_results_append(
        f"free-scalar-invert: deriv {tag} matches DFT reference = {err_deriv < check_eps}"
    )

# "half" is the default and the odd orders do not depend on even_deriv_kernel
f_half = fft_f * src
q.free_scalar_deriv_mom(f_half, [2, 0, 0, 0], even_deriv_kernel="half")
f_def = fft_f * src
q.free_scalar_deriv_mom(f_def, [2, 0, 0, 0])
err_def = float(
    np.max(np.abs(get_global_arr(geo, f_half) - get_global_arr(geo, f_def)))
)
assert err_def == 0.0, err_def
q.json_results_append(
    f"free-scalar-invert: even_deriv_kernel half is the default = {err_def == 0.0}"
)
f_odd_half = fft_f * src
q.free_scalar_deriv_mom(f_odd_half, [1, 0, 0, 0], even_deriv_kernel="half")
f_odd_cen = fft_f * src
q.free_scalar_deriv_mom(f_odd_cen, [1, 0, 0, 0], even_deriv_kernel="central")
err_odd_scheme = float(
    np.max(np.abs(get_global_arr(geo, f_odd_half) - get_global_arr(geo, f_odd_cen)))
)
assert err_odd_scheme == 0.0, err_odd_scheme
q.json_results_append(
    f"free-scalar-invert: odd order independent of even_deriv_kernel = {err_odd_scheme == 0.0}"
)
is_bad_kernel_rejected = False
try:
    f_bad = fft_f * src
    q.free_scalar_deriv_mom(f_bad, [1, 0, 0, 0], even_deriv_kernel="bad")
except Exception:
    is_bad_kernel_rejected = True
assert is_bad_kernel_rejected
q.json_results_append(
    f"free-scalar-invert: bad even_deriv_kernel is rejected = {is_bad_kernel_rejected}"
)

# self-conjugate momentum k_mu = pi: the odd factor i sin(k_mu) vanishes there,
# so the mode is dropped for odd orders and kept for even orders
for mu, tag in [(0, "d0"), (3, "d3")]:
    edge_mask = mk_self_conjugate_mask(latt_size, mu, zero_twist)
    assert bool(np.any(edge_mask)), "no self-conjugate momentum on this lattice"
    deriv_odd = [1 if i == mu else 0 for i in range(4)]
    f_odd = fft_f * src
    q.free_scalar_deriv_mom(f_odd, deriv_odd)
    err_odd_edge = float(np.max(np.abs(get_global_arr(geo, f_odd)[edge_mask])))
    assert err_odd_edge == 0.0, (mu, err_odd_edge)
    q.json_results_append(
        f"free-scalar-invert: deriv {tag} drops self-conjugate mode = {err_odd_edge == 0.0}"
    )
    deriv_even = [2 if i == mu else 0 for i in range(4)]
    f_even = fft_f * src
    q.free_scalar_deriv_mom(f_even, deriv_even)
    even_edge = get_global_arr(geo, f_even)[edge_mask]
    err_even_edge = float(np.max(np.abs(even_edge + 4.0 * src_mom_arr[edge_mask])))
    assert err_even_edge < check_eps, (mu, err_even_edge)
    q.json_results_append(
        f"free-scalar-invert: deriv {tag}^2 keeps self-conjugate mode = {err_even_edge < check_eps}"
    )

# with even_deriv_kernel="central" the even kernel is the square of the odd
# factor, so the self-conjugate mode is dropped for even orders as well
for mu, tag in [(0, "d0"), (3, "d3")]:
    edge_mask = mk_self_conjugate_mask(latt_size, mu, zero_twist)
    deriv_even = [2 if i == mu else 0 for i in range(4)]
    f_even = fft_f * src
    q.free_scalar_deriv_mom(f_even, deriv_even, even_deriv_kernel="central")
    err_even_edge_c = float(np.max(np.abs(get_global_arr(geo, f_even)[edge_mask])))
    assert err_even_edge_c == 0.0, (mu, err_even_edge_c)
    q.json_results_append(
        f"free-scalar-invert: deriv {tag}^2 central drops self-conjugate mode = {err_even_edge_c == 0.0}"
    )

# local check: the odd factor i sin(k_mu) is the central difference
# (f(x+1) - f(x-1)) / 2 in the mu direction
for mu in [0, 3]:
    deriv = [0, 0, 0, 0]
    deriv[mu] = 1
    f_mom = fft_f * src
    q.free_scalar_deriv_mom(f_mom, deriv)
    out_arr = get_global_arr(geo, fft_b * f_mom)
    cd_arr = 0.5 * (np.roll(src_arr, -1, axis=mu) - np.roll(src_arr, 1, axis=mu))
    err_cd = float(np.max(np.abs(out_arr - cd_arr)))
    assert err_cd < check_eps, (mu, err_cd)
    q.json_results_append(
        f"free-scalar-invert: deriv d{mu} is central difference = {err_cd < check_eps}"
    )

# local check: the even order deriv[mu] = 2 is unchanged and equals minus the
# mu term of the lattice laplacian
for mu in [0, 3]:
    deriv = [0, 0, 0, 0]
    deriv[mu] = 2
    f_mom = fft_f * src
    q.free_scalar_deriv_mom(f_mom, deriv)
    out_arr = get_global_arr(geo, fft_b * f_mom)
    lap_mu_arr = (
        2.0 * src_arr - np.roll(src_arr, 1, axis=mu) - np.roll(src_arr, -1, axis=mu)
    )
    err_lap = float(np.max(np.abs(out_arr + lap_mu_arr)))
    assert err_lap < check_eps, (mu, err_lap)
    q.json_results_append(
        f"free-scalar-invert: deriv d{mu}^2 is minus laplacian term = {err_lap < check_eps}"
    )

# local check: with even_deriv_kernel="central", deriv[mu] = 2 is the square of
# the central difference, (f(x+2) - 2 f(x) + f(x-2)) / 4
for mu in [0, 3]:
    deriv = [0, 0, 0, 0]
    deriv[mu] = 2
    f_mom = fft_f * src
    q.free_scalar_deriv_mom(f_mom, deriv, even_deriv_kernel="central")
    out_arr = get_global_arr(geo, fft_b * f_mom)
    cd2_arr = 0.25 * (
        np.roll(src_arr, -2, axis=mu) - 2.0 * src_arr + np.roll(src_arr, 2, axis=mu)
    )
    err_cd2 = float(np.max(np.abs(out_arr - cd2_arr)))
    assert err_cd2 < check_eps, (mu, err_cd2)
    q.json_results_append(
        f"free-scalar-invert: deriv d{mu}^2 central is central difference squared = {err_cd2 < check_eps}"
    )

# composition with the free inverse: the derivative and 1 / D(k) commute
deriv_dx = [1, 0, 0, 0]
g_deriv_inv = fft_f * src
q.free_scalar_deriv_mom(g_deriv_inv, deriv_dx)
q.free_scalar_mom_invert(g_deriv_inv, mass)
sol_deriv_inv = get_global_arr(geo, fft_b * g_deriv_inv)
g_inv_deriv = fft_f * src
q.free_scalar_mom_invert(g_inv_deriv, mass)
q.free_scalar_deriv_mom(g_inv_deriv, deriv_dx)
sol_inv_deriv = get_global_arr(geo, fft_b * g_inv_deriv)
err_comm = float(np.max(np.abs(sol_deriv_inv - sol_inv_deriv)))
assert err_comm < check_eps, err_comm
q.json_results_append(
    f"free-scalar-invert: deriv commutes with free inverse = {err_comm < check_eps}"
)
q.json_results_append(
    "free-scalar-invert: deriv free sol(xg=(1,2,3,4))",
    np.array([sol_deriv_inv[1, 2, 3, 4].real, sol_deriv_inv[1, 2, 3, 4].imag]),
    check_eps,
)

# --- free_scalar_invert_deriv: the position-space entry point.  It must
#     reproduce the explicit momentum-space composition (derivative then
#     inverse) exactly.
sol_invert_deriv = q.free_scalar_invert_deriv(src, mass, deriv=deriv_dx)
sol_invert_deriv_arr = get_global_arr(geo, sol_invert_deriv)
err_invert_deriv = float(np.max(np.abs(sol_invert_deriv_arr - sol_deriv_inv)))
# the two paths build their own FFT instances, which cuFFT need not reproduce
# bit-for-bit, so this is a tolerance check rather than an exact one
assert err_invert_deriv < check_eps, err_invert_deriv
q.json_results_append(
    f"free-scalar-invert: invert_deriv matches mom composition = {err_invert_deriv < check_eps}"
)

# free_scalar_invert_deriv forwards even_deriv_kernel: with "central" the second
# derivative of the free inverse is the free inverse of the central-difference
# squared source
deriv_d2 = [2, 0, 0, 0]
sol_d2_central = q.free_scalar_invert_deriv(
    src, mass, deriv=deriv_d2, even_deriv_kernel="central"
)
cd2_src_arr = 0.25 * (
    np.roll(src_arr, -2, axis=0) - 2.0 * src_arr + np.roll(src_arr, 2, axis=0)
)
sol_d2_ref = q.free_scalar_invert(mk_field_from_global(geo, cd2_src_arr), mass)
err_d2_central = float(
    np.max(
        np.abs(get_global_arr(geo, sol_d2_central) - get_global_arr(geo, sol_d2_ref))
    )
)
assert err_d2_central < check_eps, err_d2_central
q.json_results_append(
    f"free-scalar-invert: invert_deriv d2 central = invert(central diff squared) = {err_d2_central < check_eps}"
)

# deriv=None must reduce to free_scalar_invert
sol_invert_none = q.free_scalar_invert_deriv(src, mass)
err_invert_none = float(np.max(np.abs(get_global_arr(geo, sol_invert_none) - sol_arr)))
assert err_invert_none < check_eps, err_invert_none
q.json_results_append(
    f"free-scalar-invert: invert_deriv deriv none = free_scalar_invert = {err_invert_none < check_eps}"
)

# mode_fft=0 must agree with the default mode_fft=1
sol_invert_deriv0 = q.free_scalar_invert_deriv(src, mass, deriv=deriv_dx, mode_fft=0)
err_invert_deriv0 = float(
    np.max(np.abs(get_global_arr(geo, sol_invert_deriv0) - sol_invert_deriv_arr))
)
assert err_invert_deriv0 < check_eps, err_invert_deriv0
q.json_results_append(
    f"free-scalar-invert: invert_deriv mode_fft=0 vs mode_fft=1 = {err_invert_deriv0 < check_eps}"
)
q.json_results_append(
    "free-scalar-invert: invert_deriv sol(xg=(1,2,3,4))",
    np.array(
        [sol_invert_deriv_arr[1, 2, 3, 4].real, sol_invert_deriv_arr[1, 2, 3, 4].imag]
    ),
    check_eps,
)

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
