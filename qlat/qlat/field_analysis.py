import qlat_utils as q

import functools
import numpy as np
from .geometry import Geometry
from .c import FieldRealD, mk_fft

@functools.lru_cache(maxsize=16)
@q.timer
def mk_shift_xg_idx_arr(total_site, xg_shift):
    """
    only work without comm
    `total_site` and `xg_shift` should be tuple of ints.
    """
    total_site = np.array(total_site, dtype=np.int32)
    xg_shift = np.array(xg_shift, dtype=np.int32)
    total_volume = np.prod(total_site)
    xg_idx_arr = np.arange(total_volume).reshape(tuple(reversed(total_site))).T
    geo = Geometry(q.Coordinate(total_site))
    xg_arr = geo.xg_arr()
    xg_arr = (xg_arr + xg_shift) % total_site
    return xg_idx_arr[tuple(xg_arr.T)]

@q.timer
def smear_field_step_local(field, coef, n_steps=1):
    """
    Only work without comm
    Not very efficient
    """
    if n_steps == 0:
        return field.copy()
    geo = field.geo
    assert geo.num_node() == 1
    total_site = geo.total_site
    xg_shift_list = [
        [ 0, 0, 0, 1],
        [ 0, 0, 1, 0],
        [ 0, 1, 0, 0],
        [ 1, 0, 0, 0],
        [ 0, 0, 0, -1],
        [ 0, 0, -1, 0],
        [ 0, -1, 0, 0],
        [ -1, 0, 0, 0],
    ]
    n_dirs = len(xg_shift_list)
    xg_idx_arr_list = [ mk_shift_xg_idx_arr(tuple(total_site), tuple(xg_shift)) for xg_shift in xg_shift_list ]
    new_field = field.copy()
    coef_dir = coef / n_dirs / (1.0 - coef)
    for k in range(n_steps):
        for xg_idx_arr in xg_idx_arr_list:
            new_field[:] += coef_dir * field[xg_idx_arr]
        new_field *= 1.0 - coef
    return new_field

@q.timer
def smear_field_step(field, coef, n_steps=1):
    """
    Smear a density field
    Different from smearing the gauge field and measure the density.
    """
    if n_steps == 0:
        return field.copy()
    xg_shift_list = [
        [ 0, 0, 0, 1],
        [ 0, 0, 1, 0],
        [ 0, 1, 0, 0],
        [ 1, 0, 0, 0],
        [ 0, 0, 0, -1],
        [ 0, 0, -1, 0],
        [ 0, -1, 0, 0],
        [ -1, 0, 0, 0],
    ]
    xg_shift_list = [ q.Coordinate(xg) for xg in xg_shift_list ]
    n_dirs = len(xg_shift_list)
    coef_dir = coef / n_dirs / (1.0 - coef)
    new_field = field.copy()
    for k in range(n_steps):
        new_field *= 1 / coef_dir
        for xg_shift in xg_shift_list:
            f1 = field.shift(xg_shift)
            new_field += field.shift(xg_shift)
        new_field *= coef_dir * (1.0 - coef)
    return new_field

@functools.lru_cache(maxsize=128)
@q.timer_verbose
def mk_smear_mom_kernel(total_site, radius):
    r"""
    return f
    `radius` is the smear radius in lattice unit.
    `isinstance(f, FieldRealD)`
    `isinstance(total_site, tuple)`
    `f.total_site == q.Coordinate(total_site)`
    #
    f[:] == $G$
    #
    \ba
    G(\sigma, k) =
    \exp
    \Big(
    - \frac{\sigma^2}{2} \sum_\mu \sin^2\big(\frac{k_\mu}{2}\big)
    \Big)
    \ea
    where k[mu] = 2 pi * n[mu] / total_site[mu]
    """
    assert isinstance(total_site, tuple)
    total_site = q.Coordinate(total_site)
    geo = Geometry(total_site)
    f = FieldRealD(geo, 1)
    q.set_zero(f)
    total_site_arr = total_site.to_numpy()
    xg_arr = geo.xg_arr()
    assert f[:, 0].shape == xg_arr[:, 0].shape
    if radius == np.inf:
        sel = True
        for i in range(4):
            sel = sel & (xg_arr[:, i] == 0)
        f[sel] = 1
    else:
        gg_arr = (2 * np.pi / total_site_arr / 2) * xg_arr[:, :]
        gg_arr = np.sum(np.sin(gg_arr)**2, axis=-1)
        gg_arr = (2 / 4 * radius**2) * gg_arr
        gg_arr = np.exp(-gg_arr)
        f[:] = gg_arr[:, None]
    return f

@functools.lru_cache(maxsize=128)
@q.timer_verbose
def mk_spatial_smear_mom_kernel(total_site, radius):
    r"""
    return f
    `radius` is the smear radius in lattice unit.
    `isinstance(f, FieldRealD)`
    `isinstance(total_site, tuple)`
    `f.total_site == q.Coordinate(total_site)`
    #
    f[:] == $G$
    #
    \ba
    G(\sigma, k) =
    \exp
    \Big(
    - \frac{2\sigma^2}{3} \sum_i \sin^2\big(\frac{k_i}{2}\big)
    \Big)
    \ea
    where k[i] = 2 pi * n[i] / total_site[i]
    """
    assert isinstance(total_site, tuple)
    total_site = q.Coordinate(total_site)
    geo = Geometry(total_site)
    f = FieldRealD(geo, 1)
    q.set_zero(f)
    total_site_arr = total_site.to_numpy()
    xg_arr = geo.xg_arr()
    if radius == np.inf:
        sel = True
        for i in range(3):
            sel = sel & (xg_arr[:, i] == 0)
        f[sel] = 1
    else:
        gg_arr = (2 * np.pi / total_site_arr[:3] / 2) * xg_arr[:, :3]
        gg_arr = np.sum(np.sin(gg_arr)**2, axis=-1)
        gg_arr = (2 / 3 * radius**2) * gg_arr
        gg_arr = np.exp(-gg_arr)
        f[:] = gg_arr[:, None]
    return f

@q.timer
def smear_field(field, radius, *, is_only_spatial=False):
    r"""
    return smeared_field
    field must at least be complex type.
    #
    $$
    \ba
    f_\text{smear}(x)
    \approx
    \frac{
    \sum_y f(y) \exp( - (x - y)^2 / (2 r^2) )
    }{
    \sum_y \exp( - y^2 / (2 r^2) )
    }
    \ea
    $$
    """
    total_site = field.geo.total_site
    if is_only_spatial:
        fk = mk_spatial_smear_mom_kernel(total_site.to_tuple(), radius)
    else:
        fk = mk_smear_mom_kernel(total_site.to_tuple(), radius)
    fft1 = mk_fft(True, is_only_spatial=is_only_spatial, is_normalizing=True)
    fft2 = mk_fft(False, is_only_spatial=is_only_spatial, is_normalizing=True)
    smeared_field = field.copy()
    ftmp = fft1 * field
    ftmp *= fk
    ftmp = fft2 * ftmp
    smeared_field @= ftmp
    return smeared_field
