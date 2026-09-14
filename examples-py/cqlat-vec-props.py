#!/usr/bin/env python3

# Tests for the cqlat vector_utils interface (qlat/cqlat/vector_utils.cpp):
#   diff_gauge, load_gwu_link, save_gwu_prop, load_gwu_prop,
#   save_gwu_noiP, load_gwu_noiP, diff_prop, random_point_src,
#   make_point_prop, make_volume_src, local_sequential_source, meson_corr,
#   corr_dat_create, corr_dat_info, prop4d_conj, prop4d_src_gamma,
#   prop4d_sink_gamma, load_qlat_link, save_qlat_prop
#
# None of these functions has a Python wrapper, so they are called through the
# cqlat module: ``import qlat.c as qc`` then ``qc.<name>(...)``.  ``q.Prop``
# (a WilsonMatrix field), ``q.GaugeField``, ``q.FieldRealD`` and
# ``q.FieldComplexD`` are passed as the corresponding C++ arguments.
#
# The checks are:
#   * diff_gauge / diff_prop (which only print a comparison) are called on
#     equal copies and on a slightly perturbed copy,
#   * save_gwu_prop / load_gwu_prop are round tripped (the gwu writer only
#     stores single precision),
#   * save_gwu_noiP / load_gwu_noiP are round tripped exactly against a numpy
#     re-implementation of the saved phase,
#   * load_gwu_link reads a hand-built big-endian file (there is no exported
#     gwu link writer),
#   * random_point_src / make_point_prop / make_volume_src are checked with a
#     global gather of the local buffers,
#   * local_sequential_source is checked against a numpy time-slice mask,
#   * meson_corr / corr_dat_create / corr_dat_info are checked through the
#     files they write,
#   * prop4d_conj / prop4d_src_gamma / prop4d_sink_gamma are checked against an
#     independent numpy re-implementation of utils_corr_prop.h,
#   * save_qlat_prop (which actually writes a gauge field) / load_qlat_link are
#     round tripped.

import gc
import os

import numpy as np

import qlat as q
import qlat.c as qc

# ------------------------------------------------------------------ helpers

def local_indices(geo, latt_size):
    # local geometry x global coordinates x flat (C order) global index
    xs = np.array(
        [
            list(geo.coordinate_g_from_l(geo.coordinate_from_index(index)))
            for index in range(geo.local_volume)
        ],
        dtype=np.int64,
    )
    idx = np.ravel_multi_index((xs[:, 0], xs[:, 1], xs[:, 2], xs[:, 3]), latt_size)
    return xs, idx

def set_prop(prop, xs, f_g):
    # f_g is a global array with shape latt_size + (12, 12)
    np.asarray(prop)[:, 0] = f_g[xs[:, 0], xs[:, 1], xs[:, 2], xs[:, 3]]

def mk_prop(geo, xs, f_g):
    prop = q.Prop(geo)
    set_prop(prop, xs, f_g)
    return prop

def set_gauge(gf, xs, f_g):
    # f_g is a global array with shape latt_size + (4, 3, 3)
    np.asarray(gf)[:, :, :, :] = f_g[xs[:, 0], xs[:, 1], xs[:, 2], xs[:, 3]]

def mk_gauge(geo, xs, f_g):
    gf = q.GaugeField(geo)
    set_gauge(gf, xs, f_g)
    return gf

def prop_global(prop, idx, n_site):
    # gather the local (12, 12) matrices to the full lattice
    arr = np.asarray(prop)[:, 0]
    res = np.zeros((n_site, 12, 12), dtype=np.complex128)
    res[idx] = arr
    return q.glb_sum(res)

def ref_ga_matrices():
    # (g, ind) represents the 4x4 matrix M with M[i, ind[i]] = g[i].
    # This mirrors ga_M / set_GAM / ga_matrices_cps of utils_gammas.h.
    def mmul(x, y):
        gx, ix = x
        gy, iy = y
        return (gx * gy[ix], iy[ix])
    #
    unit = (np.ones(4, dtype=np.complex128), np.arange(4))
    i4 = np.arange(4)
    a = np.array([-1.0, 1.0, 1.0, -1.0])
    ga_i = [
        (1j * (1.0 - 2.0 * (i4 // 2)), 3 - i4),
        (a.astype(np.complex128), 3 - i4),
        (1j * (-1.0 * a), (i4 + 2) % 4),
        (np.ones(4, dtype=np.complex128), (i4 + 2) % 4),
        ((1.0 - 2.0 * (i4 // 2)).astype(np.complex128), i4),
    ]
    ga0 = [unit] + ga_i
    ga = [[mmul(ga0[i], ga0[j]) for j in range(6)] for i in range(6)]
    gL = []
    for i in range(6):
        gL.append(ga[0][i])
    for i in range(2, 6):
        gL.append(ga[1][i])
    for i in range(3, 6):
        gL.append(ga[2][i])
    for i in range(4, 6):
        gL.append(ga[3][i])
    for i in range(5, 6):
        gL.append(ga[4][i])
    return gL

def ref_src_gamma(arr, gm, conj):
    # utils_corr_prop.h prop4d_src_gammaT<dir=0>
    gg, ind = gm
    if conj:
        gg = np.conj(gg)
        src = np.conj(arr)
    else:
        src = arr
    out = arr.copy()
    for s in range(4):
        for c0 in range(3):
            for d0 in range(4):
                out[:, s * 3 + c0, ind[d0] * 3 + np.arange(3)] = (
                    gg[d0] * src[:, s * 3 + c0, d0 * 3 + np.arange(3)]
                )
    return out

def ref_sink_gamma(arr, gm, conj):
    # utils_corr_prop.h prop4d_src_gammaT<dir=1>
    gg, ind = gm
    if conj:
        gg = np.conj(gg)
        src = np.conj(arr)
    else:
        src = arr
    out = arr.copy()
    for s in range(4):
        for c0 in range(3):
            for d0 in range(4):
                out[:, d0 * 3 + c0, s * 3 + np.arange(3)] = (
                    gg[d0] * src[:, ind[d0] * 3 + c0, s * 3 + np.arange(3)]
                )
    return out

def ref_prop4d_conj(arr, rotate):
    # utils_corr_prop.h prop4d_conj, reproduced bug-for-bug: the C++ writes the
    # same destination entry for every c1 (the sink column index uses c0 instead
    # of c1), so the last c1 == 2 wins and only 48 of the 144 entries per site
    # are ever written; the nesting below reproduces that.  The header carries a
    # matching "KNOWN BUG" note and these checks pin the current behavior, so
    # they have to be updated together with any fix.
    out = arr.copy()
    for c0 in range(3):
        for d0 in range(4):
            for c1 in range(3):
                for d1 in range(4):
                    if rotate == 1:
                        out[:, d0 * 3 + c0, d1 * 3 + c0] = np.conj(
                            arr[:, d1 * 3 + c1, d0 * 3 + c0]
                        )
                    else:
                        out[:, d0 * 3 + c0, d1 * 3 + c0] = np.conj(
                            arr[:, d0 * 3 + c0, d1 * 3 + c1]
                        )
    return out

def ref_gwu_link():
    # utils_io_vec.h load_gwu_link: for each site and direction
    #   gf[x][(dir*3+c0)*3+c1] = link[((dir*3+c1)*3+c0)*2 + 0]
    #                          + i * link[((dir*3+c1)*3+c0)*2 + 1]
    # and link block v is the constant float(v).
    val = np.arange(72, dtype=np.float64)
    res = np.zeros((4, 3, 3), dtype=np.complex128)
    for d in range(4):
        for c0 in range(3):
            for c1 in range(3):
                v = ((d * 3 + c1) * 3 + c0) * 2
                res[d, c0, c1] = val[v] + 1j * val[v + 1]
    return res

# --------------------------------------------------------------- the test

q.begin_with_mpi([[2, 1, 1, 1], [1, 1, 1, 1]])

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))
q.json_results_append(f"cqlat-vec-props: geo.show()={geo.show()}")
latt_size = tuple(int(geo.total_site[i]) for i in range(4))
n_site = int(np.prod(latt_size))
nt = latt_size[3]
xs, idx_l = local_indices(geo, latt_size)
# np.ravel_multi_index uses C order, so the t coordinate is the fastest axis
t_of_site = np.arange(n_site) % nt
rs = q.RngState("cqlat-vec-props")
gL = ref_ga_matrices()

prop_g = rs.split("prop").u_rand_arr(latt_size + (12, 12)).astype(np.complex128)
prop2_g = rs.split("prop2").u_rand_arr(latt_size + (12, 12)).astype(np.complex128)
gf_g = rs.split("gauge").u_rand_arr(latt_size + (4, 3, 3)).astype(np.complex128)
seq_g = rs.split("seq").u_rand_arr(latt_size + (12, 12)).astype(np.complex128)
prop_gf = prop_g.reshape(n_site, 12, 12)
prop2_gf = prop2_g.reshape(n_site, 12, 12)
seq_gf = seq_g.reshape(n_site, 12, 12)
eye12 = np.eye(12, dtype=np.complex128)

# --- diff_gauge / diff_prop -------------------------------------------------
gf_a = mk_gauge(geo, xs, gf_g)
gf_b = mk_gauge(geo, xs, gf_g)
qc.diff_gauge(gf_a, gf_b)
assert np.array_equal(np.asarray(gf_a), np.asarray(gf_b))
gf_c = mk_gauge(geo, xs, gf_g)
np.asarray(gf_c)[0, 0, 0, 0] += 1e-9
qc.diff_gauge(gf_a, gf_c)
err = q.glb_sum(float(np.abs(np.asarray(gf_a) - np.asarray(gf_c)).max()))
assert 0.0 < err < 1e-6, err
q.json_results_append("cqlat-vec-props: diff_gauge perturbed diff", err, 1e-9)

p_a = mk_prop(geo, xs, prop_g)
p_b = mk_prop(geo, xs, prop_g)
qc.diff_prop(p_a, p_b)
assert np.array_equal(np.asarray(p_a), np.asarray(p_b))
p_c = mk_prop(geo, xs, prop_g)
np.asarray(p_c)[0, 0, 0, 0] *= 1.0 + 1e-9
qc.diff_prop(p_a, p_c)
err = q.glb_sum(float(np.abs(np.asarray(p_a) - np.asarray(p_c)).max()))
assert err > 0.0
q.json_results_append("cqlat-vec-props: diff_prop perturbed diff", err, 1e-8)

# --- save_gwu_prop / load_gwu_prop (single precision round trip) ------------
gwu_prop_path = "vec-props-tmp.gwu-prop"
qc.save_gwu_prop(p_a, gwu_prop_path)
p_gwu = q.Prop(geo)
qc.load_gwu_prop(p_gwu, gwu_prop_path)
dmax = q.glb_sum(float(np.abs(np.asarray(p_gwu) - np.asarray(p_a)).max()))
smax = q.glb_sum(float(np.abs(np.asarray(p_a)).max()))
err_gwu = dmax / smax
assert err_gwu < 1e-5, err_gwu
q.json_results_append("cqlat-vec-props: save/load_gwu_prop rel err", err_gwu, 1e-5)

# --- save_gwu_noiP / load_gwu_noiP -----------------------------------------
gwu_noi_path = "vec-props-tmp.gwu-noi"
qc.save_gwu_noiP(p_a, gwu_noi_path)
p_noi = q.Prop(geo)
qc.load_gwu_noiP(p_noi, gwu_noi_path)
arr_a = np.asarray(p_a)[:, 0]
arr_noi = np.asarray(p_noi)[:, 0]
# The saved value is the raw prop(x)(0,0) element gated by the 1-norm of the
# whole 12x12 matrix -- it is deliberately *not* normalised to a unit modulus
# phase (see the KNOWN QUIRK note on save_gwu_noiP in utils_io_vec.h), so the
# reference below reproduces that rule exactly.
exp_noi = np.zeros_like(arr_a)
n_noi = 0
for i in range(geo.local_volume):
    if float(np.abs(arr_a[i]).sum()) > 1e-8:
        n_noi += 1
        for d0 in range(12):
            exp_noi[i, d0, d0] = arr_a[i, 0, 0]
err_noi = float(np.abs(arr_noi - exp_noi).max())
assert n_noi == geo.local_volume, n_noi
assert err_noi == 0.0, err_noi
q.json_results_append("cqlat-vec-props: save/load_gwu_noiP max err", err_noi, 1e-12)

# --- load_gwu_link ----------------------------------------------------------
# There is no exported gwu link writer, so build the file with numpy.  The
# reader (read_kentucky_vector, Nvec = 4*9*2 = 72, gN = 18, dsize = 8) reads
# 72 consecutive blocks of vol = prod(total_site) = 512 big-endian float64.
gwu_link_path = "vec-props-tmp.gwu-link"
if q.get_id_node() == 0:
    np.repeat(np.arange(72, dtype=np.float64), n_site).astype(">f8").tofile(
        gwu_link_path
    )
q.sync_node()
gf_link = q.GaugeField(geo)
qc.load_gwu_link(gf_link, gwu_link_path)
arr_link = np.asarray(gf_link)
exp_link = ref_gwu_link()
assert np.array_equal(arr_link, np.broadcast_to(exp_link, arr_link.shape))
q.json_results_append(
    "cqlat-vec-props: load_gwu_link [0,0,0,0].real",
    float(arr_link[0, 0, 0, 0].real),
    1e-12,
)
q.json_results_append(
    "cqlat-vec-props: load_gwu_link [0,0,0,0].imag",
    float(arr_link[0, 0, 0, 0].imag),
    1e-12,
)
q.json_results_append(
    "cqlat-vec-props: load_gwu_link sum abs",
    q.glb_sum(float(np.abs(arr_link).sum())),
    1e-8,
)

# --- random_point_src -------------------------------------------------------
p_rnd = q.Prop(geo)
qc.random_point_src(p_rnd, 20240913)
rnd_g = prop_global(p_rnd, idx_l, n_site)
nz = int(np.count_nonzero(np.abs(rnd_g) > 1e-12))
assert nz == 12, nz
sum_abs = float(np.abs(rnd_g).sum())
assert abs(sum_abs - 12.0) < 1e-12, sum_abs
site = int(np.argmax(np.abs(rnd_g).sum(axis=(1, 2))))
assert np.array_equal(rnd_g[site], eye12)
coord = tuple(int(v) for v in np.unravel_index(site, latt_size))
q.json_results_append(
    "cqlat-vec-props: random_point_src nonzero count", float(nz), 1e-9
)
q.json_results_append("cqlat-vec-props: random_point_src sum abs", sum_abs, 1e-9)
q.json_results_append(f"cqlat-vec-props: random_point_src site = {coord}")

# --- make_point_prop --------------------------------------------------------
sp = [1, 2, 3, 4]
p_pt = q.Prop(geo)
qc.make_point_prop(p_pt, sp)
pt_g = prop_global(p_pt, idx_l, n_site)
nz = int(np.count_nonzero(np.abs(pt_g) > 1e-12))
assert nz == 12, nz
site = int(np.argmax(np.abs(pt_g).sum(axis=(1, 2))))
assert tuple(int(v) for v in np.unravel_index(site, latt_size)) == tuple(sp)
assert np.array_equal(pt_g[site], eye12)
q.json_results_append("cqlat-vec-props: make_point_prop nonzero count", float(nz), 1e-9)
q.json_results_append(
    "cqlat-vec-props: make_point_prop sum abs", float(np.abs(pt_g).sum()), 1e-9
)

# --- make_volume_src --------------------------------------------------------
# seed == -1 gives a deterministic unit diagonal on every selected site
p_vol = q.Prop(geo)
qc.make_volume_src(p_vol, -1)
vol_g = prop_global(p_vol, idx_l, n_site)
err_vol = float(np.abs(vol_g - np.broadcast_to(eye12, vol_g.shape)).max())
assert err_vol == 0.0, err_vol
q.json_results_append("cqlat-vec-props: make_volume_src(-1) max err", err_vol, 1e-12)

# tini restricts the source to one time slice
p_vol_t = q.Prop(geo)
qc.make_volume_src(p_vol_t, -1, 0, 0, 3)
vt_g = prop_global(p_vol_t, idx_l, n_site)
mask_t = t_of_site == 3
nz = int(np.count_nonzero(np.abs(vt_g) > 1e-12))
assert nz == int(np.sum(mask_t)) * 12, nz
assert np.array_equal(vt_g[mask_t], np.broadcast_to(eye12, vt_g[mask_t].shape))
assert np.count_nonzero(vt_g[~mask_t]) == 0
q.json_results_append(
    "cqlat-vec-props: make_volume_src(-1, tini=3) nonzero", float(nz), 1e-9
)

# a random source is a diagonal matrix with a unit modulus phase per site
p_vol_r = q.Prop(geo)
qc.make_volume_src(p_vol_r, 777)
vr_g = prop_global(p_vol_r, idx_l, n_site)
ph = vr_g[:, 0, 0]
assert float(np.abs(np.abs(ph) - 1.0).max()) < 1e-12
for d0 in range(1, 12):
    assert np.array_equal(vr_g[:, d0, d0], ph)
off = vr_g.copy()
for d0 in range(12):
    off[:, d0, d0] = 0.0
assert np.count_nonzero(off) == 0
nz_r = int(np.count_nonzero(np.abs(vr_g) > 1e-12))
assert nz_r == n_site * 12, nz_r
q.json_results_append(
    "cqlat-vec-props: make_volume_src(777) nonzero", float(nz_r), 1e-9
)
q.json_results_append(
    "cqlat-vec-props: make_volume_src(777) sum abs", float(np.abs(vr_g).sum()), 1e-6
)

# mix_color / mix_spin select which spin-color blocks are filled
p_vol_m1 = q.Prop(geo)
qc.make_volume_src(p_vol_m1, 777, 1, 0)
nz_m1 = int(np.count_nonzero(np.abs(prop_global(p_vol_m1, idx_l, n_site)) > 1e-12))
assert nz_m1 == n_site * 36, nz_m1
p_vol_m2 = q.Prop(geo)
qc.make_volume_src(p_vol_m2, 777, 1, 1)
nz_m2 = int(np.count_nonzero(np.abs(prop_global(p_vol_m2, idx_l, n_site)) > 1e-12))
assert nz_m2 == n_site * 72, nz_m2
q.json_results_append(
    "cqlat-vec-props: make_volume_src mix_color nonzero", float(nz_m1), 1e-9
)
q.json_results_append(
    "cqlat-vec-props: make_volume_src mix_color+spin nonzero", float(nz_m2), 1e-9
)

# --- local_sequential_source ------------------------------------------------
p_seq = mk_prop(geo, xs, seq_g)
tseq = [0, 3, 5]
mask_seq = np.isin(t_of_site, tseq)
res_seq = q.Prop(geo)
qc.local_sequential_source(res_seq, p_seq, tseq, -1)
res_g = prop_global(res_seq, idx_l, n_site)
assert np.array_equal(res_g[mask_seq], seq_gf[mask_seq])
assert np.count_nonzero(res_g[~mask_seq]) == 0
q.json_results_append("cqlat-vec-props: local_sequential_source max err", 0.0, 1e-12)

# gammai = 4 is gL[4] = ga_cps.ga[0][4]; the res is multiplied on the sink by
# that gamma (utils_corr_prop.h prop4d_sink_gamma)
res_seq4 = q.Prop(geo)
qc.local_sequential_source(res_seq4, p_seq, tseq, 4)
r4_g = prop_global(res_seq4, idx_l, n_site)
exp_seq4 = ref_sink_gamma(res_g, gL[4], False)
assert np.array_equal(r4_g, exp_seq4)
q.json_results_append(
    "cqlat-vec-props: local_sequential_source gammai=4 max err", 0.0, 1e-12
)

# --- prop4d_conj / prop4d_src_gamma / prop4d_sink_gamma ---------------------
# check that g0 == 0 is the identity (ga[0][0] == unit)
p_gam = mk_prop(geo, xs, prop2_g)
qc.prop4d_src_gamma(p_gam, 0)
got = prop_global(p_gam, idx_l, n_site)
assert np.array_equal(got, prop2_gf)
q.json_results_append("cqlat-vec-props: prop4d_src_gamma g0=0 max err", 0.0, 1e-12)

for g0 in [1, 4]:
    p_gam = mk_prop(geo, xs, prop2_g)
    qc.prop4d_src_gamma(p_gam, g0)
    got = prop_global(p_gam, idx_l, n_site)
    exp = ref_src_gamma(prop2_gf, gL[g0], False)
    err = float(np.abs(got - exp).max())
    assert err == 0.0, (g0, err)
    q.json_results_append(
        f"cqlat-vec-props: prop4d_src_gamma g0={g0} max err", err, 1e-12
    )

p_gam = mk_prop(geo, xs, prop2_g)
qc.prop4d_src_gamma(p_gam, 1, 1)
got = prop_global(p_gam, idx_l, n_site)
err = float(np.abs(got - ref_src_gamma(prop2_gf, gL[1], True)).max())
assert err == 0.0, err
q.json_results_append("cqlat-vec-props: prop4d_src_gamma g0=1 conj max err", err, 1e-12)

for g0 in [1, 4]:
    for conj in [0, 1]:
        p_gam = mk_prop(geo, xs, prop2_g)
        qc.prop4d_sink_gamma(p_gam, g0, conj)
        got = prop_global(p_gam, idx_l, n_site)
        exp = ref_sink_gamma(prop2_gf, gL[g0], bool(conj))
        err = float(np.abs(got - exp).max())
        assert err == 0.0, (g0, conj, err)
        q.json_results_append(
            f"cqlat-vec-props: prop4d_sink_gamma g0={g0} conj={conj} max err",
            err,
            1e-12,
        )

for rotate in [0, 1]:
    p_gam = mk_prop(geo, xs, prop2_g)
    qc.prop4d_conj(p_gam, rotate)
    got = prop_global(p_gam, idx_l, n_site)
    exp = ref_prop4d_conj(prop2_gf, rotate)
    err = float(np.abs(got - exp).max())
    assert err == 0.0, (rotate, err)
    q.json_results_append(
        f"cqlat-vec-props: prop4d_conj rotate={rotate} max err", err, 1e-12
    )

# --- meson_corr -------------------------------------------------------------
meson_path = "vec-props-tmp-meson.corr"
qc.meson_corr(
    p_a, p_seq, meson_path, 0, 0, 0, 1, "cqlat-vec-props meson info", 1, [0, 0, 0, 0]
)
q.sync_node()
size_meson = q.glb_sum(
    float(os.path.getsize(meson_path)) if q.get_id_node() == 0 else 0.0
)
assert size_meson > 0.0
q.json_results_append("cqlat-vec-props: meson_corr file size", size_meson, 1e-9)

# --- corr_dat_create / corr_dat_info ----------------------------------------
corr_path = "vec-props-tmp-corr.dat"
info1 = "cqlat-vec-props corr info one"
info2 = "cqlat-vec-props corr info two"
qc.corr_dat_create(corr_path, "8 2", "t p", info1)
q.sync_node()
size1 = q.glb_sum(float(os.path.getsize(corr_path)) if q.get_id_node() == 0 else 0.0)
assert size1 > 0.0
if q.get_id_node() == 0:
    with open(corr_path, "rb") as fp:
        raw = fp.read()
    assert b"BEGIN_Corr_HEAD" in raw
    assert info1.encode() in raw
qc.corr_dat_info(corr_path, info2)
q.sync_node()
size2 = q.glb_sum(float(os.path.getsize(corr_path)) if q.get_id_node() == 0 else 0.0)
if q.get_id_node() == 0:
    with open(corr_path, "rb") as fp:
        raw = fp.read()
    assert info1.encode() in raw
    assert info2.encode() in raw
assert size2 >= size1
q.json_results_append("cqlat-vec-props: corr_dat_create size", size1, 1e-9)
q.json_results_append("cqlat-vec-props: corr_dat_info size", size2, 1e-9)

# --- save_qlat_prop / load_qlat_link ----------------------------------------
# save_qlat_prop is a misnomer: it forwards to save_qlat_link and writes a
# gauge field (double precision, multiplicity 4), not a propagator.  The
# matching propagate writes/reads (C++ save_qlat_prop / load_qlat_prop for
# Propagator4d) are not exported, so a Propagator cannot be written with the
# qlat format from Python; the cqlat export carries the same note.
qlat_path = "vec-props-tmp.qlat-link"
gf_q = mk_gauge(geo, xs, gf_g)
qc.save_qlat_prop(gf_q, qlat_path)
q.sync_node()
gf_q2 = q.GaugeField(geo)
qc.load_qlat_link(gf_q2, qlat_path)
err_ql = float(np.abs(np.asarray(gf_q) - np.asarray(gf_q2)).max())
assert err_ql == 0.0, err_ql
q.json_results_append(
    "cqlat-vec-props: save_qlat_prop/load_qlat_link max err", err_ql, 1e-12
)
q.json_results_append(
    "cqlat-vec-props: load_qlat_link sum abs",
    q.glb_sum(float(np.abs(np.asarray(gf_q2)).sum())),
    1e-8,
)

# --- cleanup ----------------------------------------------------------------
del p_gwu, p_noi, p_rnd, p_pt, p_vol, p_vol_t, p_vol_r, p_vol_m1, p_vol_m2
del p_seq, res_seq, res_seq4, p_gam, gf_a, gf_b, gf_c, gf_link, gf_q, gf_q2
gc.collect()
for fn in [
    gwu_prop_path,
    gwu_noi_path,
    gwu_link_path,
    meson_path,
    corr_path,
    qlat_path,
]:
    if q.get_id_node() == 0:
        if os.path.exists(fn):
            os.remove(fn)
q.sync_node()

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
