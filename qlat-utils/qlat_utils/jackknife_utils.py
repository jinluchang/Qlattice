"""
Module ``qlat_utils.jackknife_utils``
=====================================\n
Jackknife resampling and error estimation: the plain jackknife, the
super-jackknife, the randomized jackknife-bootstrap hybrid, and the unified
``g_*`` interface with its shared ``default_g_jk_kwargs``.  Each family also has
collective MPI variants, which distribute the input and/or the output between
the nodes, implemented with ``mpi4py``.\n
The generic helpers (the type tuples, ``use_kwargs``, ``average``,
``block_data``, ``filter_np_results``, ``fsqr``, ``fsqrt``, ``qnorm`` and
``NewDictValues``) are imported from ``qlat_utils.data``.\n
Documentation: ``docs/qlat-utils/qlat_jackknife_utils.md``\n
.. note:: Update the documentation when updating this source file.
"""

import math
import numpy as np

class q:
    from qlat_utils.data import (
        int_types,
        real_types,
        number_types,
        use_kwargs,
        average,
        block_data,
        filter_np_results,
        fsqr,
        fsqrt,
        qnorm,
        NewDictValues,
    )
    from qlat_utils.timer import (
        timer,
        displayln_info,
    )
    from qlat_utils.rng_state import (
        RngState,
    )

# ----------

def jackknife(data_list, *, eps=1):
    r"""
    Return jk[i] = avg - \frac{eps}{N} (v[i] - avg)
    normal jackknife uses eps=1, scale the fluctuation by eps
    """
    is_np_arr = isinstance(data_list, np.ndarray)
    data_list_real = [d for d in data_list if d is not None]
    n = len(data_list_real)
    fac = eps / n
    avg = q.average(data_list_real)
    jks = [
        avg,
    ]
    for data in data_list:
        if data is None:
            jks.append(avg)
        else:
            jks.append(avg - fac * (data - avg))
    if is_np_arr:
        jks = np.array(jks, dtype=data_list.dtype)
    return jks

# ----

def jk_avg(jk_arr):
    val = jk_arr[0]
    return q.filter_np_results(val)

def jk_err(jk_arr, *, eps=1, block_size=1):
    r"""
    Return
    $$
    \frac{1}{eps} \sqrt{ N/(N-block_size) \sum_{i=1}^N (jk[i] - jk_avg)^2 }.
    $$
    when ``block_size=1``.
    Note: ``len(jk_arr) = N + 1``.
    Same ``eps`` as the ``eps`` used in the ``jackknife`` function.
    Does not properly honor the $(N-1)$ formula in error calculation
    if there were missing data in the original ``data_list`` in the ``jackknife`` function.
    """
    assert block_size >= 1
    avg = jk_avg(jk_arr)
    n = len(jk_arr) - 1
    if n <= 1:
        fac = 1 / abs(eps)
        val = fac * avg
        val = q.filter_np_results(val)
        return val
    if n < 2 * block_size:
        block_size = 1
    assert n > block_size
    blocks = q.block_data(jk_arr[1:], block_size)
    diff_sqr = q.average([q.fsqr(jk - avg) for jk in blocks])
    fac = math.sqrt(block_size / (n - block_size)) * n / abs(eps)
    val = fac * q.fsqrt(diff_sqr)
    val = q.filter_np_results(val)
    return val

def jk_avg_err(jk_arr, *, eps=1, block_size=1):
    return jk_avg(jk_arr), jk_err(jk_arr, eps=eps, block_size=block_size)

@q.timer
def sjackknife(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    is_hash_jk_idx=True,
    jk_idx_hash_size=None,
    rng_state=None,
    all_jk_idx=None,
    get_all_jk_idx=None,
    jk_blocking_func=None,
    eps=1,
):
    """
    Super jackknife.
    Return ``jk_arr``.
    ``len(jk_idx_list) == len(data_list)``
    ``len(jk_arr) == len(all_jk_idx)``
    ``jk_idx_list`` (after processed by ``jk_blocking_func``) should be contained in ``all_jk_idx``,
    otherwise, if ``is_hash_jk_idx`` is true, then, a hash of ``jk_idx`` will be used instead.
    Ideally, ``all_jk_idx`` should only contain distinct indices.
    However, if there are repeatations, indices appear later take precedence.
    if ``all_jk_idx`` and ``get_all_jk_idx`` are both ``None``, then a trivial
    ``all_jk_idx`` will be created based on ``jk_idx_hash_size``.
    """
    if jk_idx_hash_size is None:
        jk_idx_hash_size = 1024
    if rng_state is None:
        rng_state = q.RngState("rejk")
    rs = rng_state
    assert len(jk_idx_list) == len(data_list)
    if isinstance(data_list, np.ndarray):
        dtype = data_list.dtype
    else:
        dtype = None
    data_list_real = [d for d in data_list if d is not None]
    data_arr = np.array(data_list_real, dtype=dtype)
    if avg is None:
        avg = q.average(data_arr)
    dtype = data_arr.dtype
    jk_idx_list = [jk_idx for jk_idx, d in zip(jk_idx_list, data_list) if d is not None]
    if jk_blocking_func is not None:
        jk_idx_list = [jk_blocking_func(0, jk_idx) for jk_idx in jk_idx_list]
    n = len(data_arr)
    assert n == len(jk_idx_list)
    if all_jk_idx is None:
        if get_all_jk_idx is None:
            assert is_hash_jk_idx
            all_jk_idx = [
                "avg",
            ] + list(range(jk_idx_hash_size))
        else:
            all_jk_idx = get_all_jk_idx()
    assert all_jk_idx[0] == "avg"
    n_super_sample = len(all_jk_idx) - 1
    i_dict = dict()
    for i, jk_idx in enumerate(all_jk_idx):
        jk_idx_str = str(jk_idx)
        i_dict[jk_idx_str] = i
    i_arr = np.zeros(n, dtype=np.int32)
    for j in range(n):
        jk_idx = jk_idx_list[j]
        jk_idx_str = str(jk_idx)
        if jk_idx_str in i_dict:
            i = i_dict[jk_idx_str]
        else:
            assert is_hash_jk_idx
            rsi = rs.split(jk_idx_str)
            i = 1 + int(rsi.rand_gen() % n_super_sample)
        assert i > 0
        i_arr[j] = i
    count_dict = dict()
    for j in range(n):
        i = i_arr[j]
        if i in count_dict:
            count_dict[i] += 1
        else:
            count_dict[i] = 1
    jk_arr = np.empty(
        (
            1 + n_super_sample,
            *data_arr[0].shape,
        ),
        dtype=dtype,
    )
    jk_arr[:] = avg
    data_diff = data_arr - avg
    for j in range(n):
        i = i_arr[j]
        assert i > 0
        assert i in count_dict
        if n > count_dict[i]:
            n_b = n - count_dict[i]
            fac = -eps * np.sqrt(1 / (n * n_b))
            jk_arr[i] += fac * data_diff[j]
    return jk_arr

@q.timer
def sjackknife_distributed(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    is_hash_jk_idx=True,
    jk_idx_hash_size=None,
    rng_state=None,
    all_jk_idx=None,
    get_all_jk_idx=None,
    jk_blocking_func=None,
    eps=1,
):
    r"""
    ``sjackknife`` with distributed input and distributed output.\n
    This is a collective MPI operation: every node must call it with the same
    parameters and with its own disjoint part of the data set (``data_list``
    and ``jk_idx_list`` are the local parts); every node returns its own part
    of the ``jk_arr`` of ``sjackknife``.  Node ``r`` out of ``num_node`` nodes
    holds the samples in ``range(*get_distributed_range(len(all_jk_idx), r,
    num_node))``, so that::
        jk_arr = np.concatenate(q.get_comm().allgather(jk_local))
    is the complete data set, in the same order as the one obtained with
    ``sjackknife``.  No node needs to hold the whole data set or the whole
    result; only the small ``jk_idx`` metadata is gathered on every node.\n
    The result agrees with the one of ``sjackknife`` up to the floating-point
    roundoff, but not bit-for-bit: the average and the sums over the data set
    are reduced across the nodes, which changes the order of the
    floating-point additions.\n
    ``sjackknife_sync_node`` performs the corresponding collective operation
    for the case where every node has the whole input; ``g_mk_jk_distributed``
    dispatches to this function.
    """
    fname = "sjackknife_distributed"
    comm, id_node, num_node = get_collective_comm(fname)
    if jk_idx_hash_size is None:
        jk_idx_hash_size = 1024
    if rng_state is None:
        rng_state = q.RngState("rejk")
    (
        data_arr,
        jk_idx_list_local,
        jk_idx_list_glb,
        n,
        elem_shape,
        dtype,
    ) = get_distributed_jk_input(fname, comm, data_list, jk_idx_list)
    avg = get_distributed_avg(comm, data_arr, n, avg)
    if all_jk_idx is None:
        if get_all_jk_idx is None:
            assert is_hash_jk_idx
            all_jk_idx = [
                "avg",
            ] + list(range(jk_idx_hash_size))
        else:
            all_jk_idx = get_all_jk_idx()
    assert all_jk_idx[0] == "avg"
    total_size = len(all_jk_idx)
    i_start, i_end = get_distributed_range(total_size, id_node, num_node)
    # ``partial_arr`` holds the contribution of the local data to every sample
    partial_arr = np.zeros((total_size,) + elem_shape, dtype=dtype)
    if i_start == 0 and i_end > i_start:
        # the sample 0 is the average and must be counted only once
        partial_arr[0] = avg
    #
    n_super_sample = total_size - 1
    if jk_blocking_func is None:
        b_jk_idx_list_local = jk_idx_list_local
        b_jk_idx_list_glb = jk_idx_list_glb
    else:
        b_jk_idx_list_local = [
            jk_blocking_func(0, jk_idx) for jk_idx in jk_idx_list_local
        ]
        b_jk_idx_list_glb = [jk_blocking_func(0, jk_idx) for jk_idx in jk_idx_list_glb]
    i_dict = dict()
    for i, jk_idx in enumerate(all_jk_idx):
        i_dict[str(jk_idx)] = i
    rs = rng_state
    #
    def get_i(jk_idx):
        jk_idx_str = str(jk_idx)
        if jk_idx_str in i_dict:
            i = i_dict[jk_idx_str]
        else:
            assert is_hash_jk_idx
            i = 1 + int(rs.split(jk_idx_str).rand_gen() % n_super_sample)
        assert i > 0
        return i
    #
    count_dict = dict()
    for jk_idx in b_jk_idx_list_glb:
        i = get_i(jk_idx)
        count_dict[i] = count_dict.get(i, 0) + 1
    data_diff = data_arr - avg
    for j in range(len(data_arr)):
        i = get_i(b_jk_idx_list_local[j])
        count = count_dict[i]
        if n > count:
            fac = -eps * np.sqrt(1 / (n * (n - count)))
            partial_arr[i] += fac * data_diff[j]
    return get_reduce_scattered_jk_arr(partial_arr, comm, id_node, num_node, avg)

@q.timer
def sjackknife_sync_node(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    is_hash_jk_idx=True,
    jk_idx_hash_size=None,
    rng_state=None,
    all_jk_idx=None,
    get_all_jk_idx=None,
    jk_blocking_func=None,
    eps=1,
):
    r"""
    ``sjackknife`` as a collective MPI operation where every node has the same
    (whole) input and obtains the same complete result.\n
    This implements the ``is_sync_node=True`` option of ``sjackknife``: the
    input is split between the nodes in the order of the nodes and
    ``sjackknife_distributed`` is called on the local parts; the local results
    are then gathered with ``Allgatherv``, so that every node obtains the
    complete ``jk_arr`` of ``sjackknife``.  The result agrees with the one of
    ``sjackknife`` up to the floating-point roundoff, but not bit-for-bit.\n
    ``g_mk_jk_sync_node`` dispatches to this function.
    """
    comm, id_node, num_node = get_collective_comm("sjackknife_sync_node")
    i_start, i_end = get_distributed_range(len(data_list), id_node, num_node)
    jk_local = sjackknife_distributed(
        data_list[i_start:i_end],
        jk_idx_list[i_start:i_end],
        avg=avg,
        is_hash_jk_idx=is_hash_jk_idx,
        jk_idx_hash_size=jk_idx_hash_size,
        rng_state=rng_state,
        all_jk_idx=all_jk_idx,
        get_all_jk_idx=get_all_jk_idx,
        jk_blocking_func=jk_blocking_func,
        eps=eps,
    )
    return get_gathered_jk_arr(jk_local, comm, num_node)

@q.timer
def sjk_mk_jk_val(
    rs_tag,
    val,
    err,
    *,
    is_hash_jk_idx=True,
    jk_idx_hash_size=None,
    rng_state=None,
    all_jk_idx=None,
    get_all_jk_idx=None,
    eps=1,
):
    """
    return jk_arr
    n = n_rand_sample
    len(jk_arr) == 1 + n
    jk_arr[i] = val + err * r[i] for i in 1..n
    where r[i] ~ N(0, 1)
    """
    if jk_idx_hash_size is None:
        jk_idx_hash_size = 1024
    assert jk_idx_hash_size >= 0
    if rng_state is None:
        rng_state = q.RngState("rejk")
    rs = rng_state
    if all_jk_idx is None:
        if get_all_jk_idx is None:
            assert is_hash_jk_idx
            all_jk_idx = [
                "avg",
            ] + list(range(jk_idx_hash_size))
        else:
            all_jk_idx = get_all_jk_idx()
    assert all_jk_idx[0] == "avg"
    n_super_sample = len(all_jk_idx) - 1
    assert n_super_sample >= 0
    assert isinstance(rng_state, q.RngState)
    assert isinstance(val, q.real_types)
    assert isinstance(err, q.real_types)
    rs = rng_state.split(str(rs_tag))
    jk_arr = np.zeros((n_super_sample + 1,), dtype=np.float64)
    jk_arr[0] = val
    r_arr = rs.g_rand_arr((n_super_sample,))
    r_arr_qnorm = q.qnorm(r_arr)
    r_arr = r_arr * np.sqrt(1 / r_arr_qnorm)
    assert abs(q.qnorm(r_arr) - 1) < 1e-8
    jk_arr[1:] = val + eps * r_arr * err
    return jk_arr

def sjk_avg(jk_arr):
    return jk_avg(jk_arr)

def sjk_err(jk_arr, *, eps=1):
    r"""
    Return
    $$
    \frac{1}{eps} \sqrt{ \sum_{i=1}^N (jk[i] - jk_avg)^2 }.
    $$
    Note: ``len(jk_arr) = N + 1``.
    Same ``eps`` as the ``eps`` used in the ``jackknife`` function.
    """
    avg = jk_avg(jk_arr)
    n = len(jk_arr) - 1
    if n <= 1:
        fac = 1 / abs(eps)
        val = fac * avg
        val = q.filter_np_results(val)
        return val
    diff_sqr = q.average([q.fsqr(jk - avg) for jk in jk_arr[1:]])
    fac = math.sqrt(n) / abs(eps)
    val = fac * q.fsqrt(diff_sqr)
    val = q.filter_np_results(val)
    return val

def sjk_avg_err(jk_arr, *, eps=1):
    return sjk_avg(jk_arr), sjk_err(jk_arr, eps=eps)

# ----------

@q.timer
def mk_r_i_j_mat(
    n_rand_sample,
    jk_idx_list,
    rng_state,
    *,
    jk_blocking_func,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_use_old_rand_alg,
    i_range=None,
    jk_idx_list_for_count=None,
):
    """
    Return ``(r_arr, b_arr)``.\n
    ``r_arr`` and ``b_arr`` have shape ``(n_sample, n)`` with
    ``n_sample = i_end - i_start``, where ``i_range = (i_start, i_end,)``
    (default ``None`` means the full range ``(0, n_rand_sample,)``).
    Only the rows in ``i_range`` are computed; the per-block random streams
    do not depend on ``i_range``, so every returned row is identical to the
    corresponding row of a full-range call.\n
    ``jk_idx_list_for_count`` (default ``None`` means ``jk_idx_list``) is the
    list over which the block sizes ``b_arr`` are counted; the columns of
    ``r_arr`` and ``b_arr`` always correspond to ``jk_idx_list``.  This is used
    by ``g_mk_jk_distributed`` where the columns are the locally held
    ``jk_idx`` while the block sizes must be counted over the whole data set.
    """
    assert n_rand_sample >= 0
    if i_range is None:
        i_start, i_end = 0, n_rand_sample
    else:
        i_start, i_end = i_range
        assert 0 <= i_start <= i_end <= n_rand_sample
    n_sample = i_end - i_start
    rs = rng_state
    n = len(jk_idx_list)
    r_arr = np.empty(
        (
            n_sample,
            n,
        ),
        dtype=np.float64,
    )
    b_arr = np.empty(
        (
            n_sample,
            n,
        ),
        dtype=np.int32,
    )
    jk_idx_str_arr = np.empty(
        (
            n_sample,
            n,
        ),
        dtype=object,
    )
    jk_idx_str_set = set()
    if jk_blocking_func is None:
        is_apply_rand_sample_jk_idx_blocking_shift = False
    # ``b_arr`` counts the blocks of ``jk_idx_list_for_count`` (default
    # ``jk_idx_list``); the columns of ``r_arr`` and ``b_arr`` always
    # correspond to ``jk_idx_list``.
    is_separate_count = jk_idx_list_for_count is not None
    #
    @q.timer
    def set_jk_idx():
        if is_apply_rand_sample_jk_idx_blocking_shift:
            for i_local in range(n_sample):
                i = i_start + i_local
                for j in range(n):
                    jk_idx = jk_blocking_func(i + 1, jk_idx_list[j])
                    jk_idx_str = str(jk_idx)
                    jk_idx_str_arr[i_local, j] = jk_idx_str
                    jk_idx_str_set.add(jk_idx_str)
                count_dict = dict()
                if is_separate_count:
                    for jk_idx in jk_idx_list_for_count:
                        jk_idx_str = str(jk_blocking_func(i + 1, jk_idx))
                        count_dict[jk_idx_str] = count_dict.get(jk_idx_str, 0) + 1
                else:
                    for j in range(n):
                        jk_idx_str = jk_idx_str_arr[i_local, j]
                        count_dict[jk_idx_str] = count_dict.get(jk_idx_str, 0) + 1
                for j in range(n):
                    jk_idx_str = jk_idx_str_arr[i_local, j]
                    b_arr[i_local, j] = count_dict[jk_idx_str]
        else:
            count_dict = dict()
            for j in range(n):
                jk_idx = jk_idx_list[j]
                if jk_blocking_func is not None:
                    jk_idx = jk_blocking_func(0, jk_idx)
                jk_idx_str = str(jk_idx)
                jk_idx_str_arr[:, j] = jk_idx_str
                jk_idx_str_set.add(jk_idx_str)
                if not is_separate_count:
                    count_dict[jk_idx_str] = count_dict.get(jk_idx_str, 0) + 1
            if is_separate_count:
                for jk_idx in jk_idx_list_for_count:
                    if jk_blocking_func is not None:
                        jk_idx = jk_blocking_func(0, jk_idx)
                    jk_idx_str = str(jk_idx)
                    count_dict[jk_idx_str] = count_dict.get(jk_idx_str, 0) + 1
            if n_sample > 0:
                for j in range(n):
                    jk_idx_str = jk_idx_str_arr[0, j]
                    b_arr[:, j] = count_dict[jk_idx_str]
    #
    set_jk_idx()
    if is_use_old_rand_alg == "v1":
        assert not is_normalizing_rand_sample
        for i_local in range(n_sample):
            i = i_start + i_local
            rsi = rs.split(str(i))
            r = [
                rsi.split(jk_idx_str).g_rand_gen()
                for jk_idx_str in jk_idx_str_arr[i_local]
            ]
            for j in range(n):
                r_arr[i_local, j] = r[j]
        return r_arr, b_arr
    assert not is_use_old_rand_alg
    r_arr_dict = dict()
    #
    @q.timer
    def set_r():
        for jk_idx_str in jk_idx_str_set:
            rsi = rs.split(jk_idx_str)
            garr = rsi.g_rand_arr(n_rand_sample)
            if is_normalizing_rand_sample:
                # garr_qnorm \approx n_rand_sample
                garr_qnorm = q.qnorm(garr)
                garr = garr * np.sqrt(n_rand_sample / garr_qnorm)
                assert abs(q.qnorm(garr) / n_rand_sample - 1) < 1e-8
            r_arr_dict[jk_idx_str] = garr
    #
    set_r()
    #
    @q.timer
    def set_r_arr():
        if is_apply_rand_sample_jk_idx_blocking_shift:
            for i_local in range(n_sample):
                i = i_start + i_local
                for j in range(n):
                    jk_idx_str = jk_idx_str_arr[i_local, j]
                    garr = r_arr_dict[jk_idx_str]
                    r_arr[i_local, j] = garr[i]
        else:
            if n_sample > 0:
                for j in range(n):
                    jk_idx_str = jk_idx_str_arr[0, j]
                    garr = r_arr_dict[jk_idx_str]
                    r_arr[:, j] = garr[i_start:i_end]
    #
    set_r_arr()
    return r_arr, b_arr

def get_distributed_range(total_size, id_node, num_node):
    """
    Return ``(i_start, i_end)``, the range of the ``total_size`` indices owned
    by node ``id_node`` out of ``num_node`` nodes.\n
    The ranges are contiguous, ordered by ``id_node``, cover
    ``range(total_size)`` exactly once and differ in size by at most one::
        i_start = (total_size * id_node) // num_node
        i_end = (total_size * (id_node + 1)) // num_node
    Concatenating the parts in the order of the nodes reproduces the whole
    range, so ``np.concatenate(comm.allgather(x_local))`` is the complete
    array when ``x_local = x[i_start:i_end]``.
    """
    assert isinstance(total_size, q.int_types)
    assert 0 <= total_size
    assert isinstance(id_node, q.int_types)
    assert isinstance(num_node, q.int_types)
    assert 0 <= id_node < num_node
    i_start = (total_size * id_node) // num_node
    i_end = (total_size * (id_node + 1)) // num_node
    return i_start, i_end

def is_mpi_dtype(dtype):
    """
    Return whether ``dtype`` is a numpy dtype that MPI can use to sum buffers,
    i.e. a dtype supported by ``MPI.SUM``, like the ones used by
    ``g_mk_jk_distributed``.
    """
    dtype = np.dtype(dtype)
    if dtype.kind == "f":
        return dtype.itemsize in (4, 8)
    if dtype.kind == "c":
        return dtype.itemsize in (8, 16)
    if dtype.kind in ("i", "u"):
        return dtype.itemsize in (1, 2, 4, 8)
    return False

def get_collective_comm(tag):
    """
    Return ``(comm, id_node, num_node)`` for a collective MPI operation.\n
    ``tag`` names the caller for the error messages, e.g.
    ``"g_mk_jk_sync_node"``.\n
    ``qlat`` (used for ``q.get_comm()``) is imported only here.  The
    communicator must be initialized on the whole MPI communicator, i.e. with
    ``q.begin_with_mpi()``, ``q.begin_with_gpt()`` or ``q.begin_with_grid()``
    (or set with ``q.set_comm(...)``), so that
    ``comm.rank == q.get_id_node()`` and ``comm.size == q.get_num_node()``.
    """
    import qlat
    #
    comm = qlat.get_comm()
    if comm is None:
        raise Exception(
            f"{tag} requires the qlat communicator;"
            " use q.begin_with_mpi(), q.begin_with_gpt() or"
            " q.begin_with_grid() (or set it with q.set_comm(...))"
        )
    if comm.size != qlat.get_num_node():
        raise Exception(
            f"{tag} requires qlat to be initialized"
            f" on the whole MPI communicator, but comm.size={comm.size} and"
            f" q.get_num_node()={qlat.get_num_node()}"
        )
    assert comm.rank == qlat.get_id_node()
    return comm, comm.rank, comm.size

def get_distributed_jk_input(fname, comm, data_list, jk_idx_list):
    r"""
    Return ``(data_arr, jk_idx_list_local, jk_idx_list_glb, n, elem_shape, dtype)``
    for the distributed jackknife functions.\n
    ``data_list`` and ``jk_idx_list`` are the local parts of the data set and
    ``comm`` is the communicator from ``get_collective_comm``.  ``data_arr`` is
    the local data, with shape ``(n_local, *elem_shape)`` (``n_local == 0``
    when the node holds no data); ``jk_idx_list_local`` are the local
    ``jk_idx`` and ``jk_idx_list_glb`` is the whole list of ``jk_idx`` gathered
    on every node; ``n`` is the total number of data points; ``elem_shape`` and
    ``dtype`` describe a single data point.\n
    ``fname`` names the caller for the error messages.  The data must have a
    dtype supported by MPI, so that the results can be summed and gathered.
    """
    from mpi4py import MPI
    #
    assert len(data_list) == len(jk_idx_list)
    if isinstance(data_list, np.ndarray):
        dtype = data_list.dtype
    else:
        dtype = None
    data_list_local = [d for d in data_list if d is not None]
    jk_idx_list_local = [
        jk_idx for jk_idx, d in zip(jk_idx_list, data_list) if d is not None
    ]
    data_arr = np.array(data_list_local, dtype=dtype)
    n_local = len(data_arr)
    # the dtype and the element shape of the data; the lowest rank with data
    # provides them, so that a node without any data still knows them
    if n_local > 0:
        elem_info = (data_arr.dtype.str, tuple(data_arr.shape[1:]))
    else:
        elem_info = None
    elem_info = next(
        (info for info in comm.allgather(elem_info) if info is not None), None
    )
    if elem_info is None:
        raise Exception(f"{fname}: the distributed data set is empty")
    dtype = np.dtype(elem_info[0])
    elem_shape = tuple(elem_info[1])
    if not is_mpi_dtype(dtype):
        raise Exception(
            f"{fname}: the data must have a dtype supported by MPI so that the"
            f" results can be summed, but dtype={dtype}"
        )
    if n_local == 0:
        data_arr = np.zeros((0,) + elem_shape, dtype=dtype)
    else:
        data_arr = np.asarray(data_arr, dtype=dtype)
    n_arr = np.array([n_local], dtype=np.int64)
    comm.Allreduce(MPI.IN_PLACE, n_arr, op=MPI.SUM)
    n = int(n_arr[0])
    # only the small ``jk_idx`` metadata is gathered on every node
    jk_idx_list_glb = [
        jk_idx
        for jk_idx_list_node in comm.allgather(jk_idx_list_local)
        for jk_idx in jk_idx_list_node
    ]
    return data_arr, jk_idx_list_local, jk_idx_list_glb, n, elem_shape, dtype

def get_distributed_avg(comm, data_arr, n, avg=None):
    r"""
    Return the average of the whole distributed data set.\n
    When ``avg`` is not ``None``, it is returned as it is; otherwise the local
    sums of ``data_arr`` are summed over the nodes with ``Allreduce`` and
    divided by ``n``.  This is the mean over the whole data set, which is
    needed by the (randomized) Super-Jackknife, and it agrees with the one of
    the sequential functions up to the floating-point roundoff.
    """
    from mpi4py import MPI
    #
    if avg is not None:
        return avg
    local_sum = np.asarray(np.sum(data_arr, axis=0))
    glb_sum = np.zeros_like(local_sum)
    comm.Allreduce(local_sum, glb_sum, op=MPI.SUM)
    return q.filter_np_results(glb_sum / n)

def get_reduce_scattered_jk_arr(partial_arr, comm, id_node, num_node, avg):
    """
    Return the local part of the (randomized) Super-Jackknife data set.\n
    ``partial_arr`` has shape ``(total_size, *elem_shape)`` and holds the
    contribution of the local data to every sample, with the sample 0 equal to
    ``avg`` on the node which owns it and 0 elsewhere.  The contributions of
    the nodes are summed with ``Reduce_scatter``, so that every node obtains
    its own samples; ``avg`` is then added to the samples, since the samples
    are ``avg + sum_j (...)`` while the sample 0 is ``avg`` itself.\n
    ``partial_arr`` may be any array, in particular a non-contiguous view such
    as a column of a 2-D array; a contiguous copy is made when needed.
    """
    from mpi4py import MPI
    #
    # The buffers of the collectives must be contiguous; ``reshape(-1)`` is
    # only a view when the array is already 1-D, so copy a strided array here.
    partial_arr = np.ascontiguousarray(partial_arr)
    total_size = partial_arr.shape[0]
    elem_shape = partial_arr.shape[1:]
    elem_size = 1
    for x in elem_shape:
        elem_size *= x
    recvcounts = []
    for r in range(num_node):
        r_start, r_end = get_distributed_range(total_size, r, num_node)
        recvcounts.append((r_end - r_start) * elem_size)
    i_start, i_end = get_distributed_range(total_size, id_node, num_node)
    jk_arr = np.empty((i_end - i_start,) + elem_shape, dtype=partial_arr.dtype)
    comm.Reduce_scatter(
        partial_arr.reshape(-1),
        jk_arr.reshape(-1),
        recvcounts,
        op=MPI.SUM,
    )
    if i_end - i_start > 0 and i_end > 1:
        jk_arr[max(1 - i_start, 0) :] += avg
    return jk_arr

def get_gathered_jk_arr(jk_local, comm, num_node):
    """
    Return the complete (randomized) Super-Jackknife data set.\n
    ``jk_local`` is the local part of the result of a distributed jackknife
    function; the parts are gathered with ``Allgatherv``, in the order of the
    nodes, so that every node obtains the complete ``jk_arr``.\n
    ``jk_local`` may be any array, in particular a non-contiguous view such as
    a column of a 2-D array; a contiguous copy is made when needed.
    """
    from mpi4py import MPI
    #
    # The send buffer of the collective must be contiguous; ``reshape(-1)`` is
    # only a view when the array is already 1-D, so copy a strided array here.
    jk_local = np.ascontiguousarray(jk_local)
    n_arr = np.array([len(jk_local)], dtype=np.int64)
    comm.Allreduce(MPI.IN_PLACE, n_arr, op=MPI.SUM)
    total_size = int(n_arr[0])
    elem_size = 1
    for x in jk_local.shape[1:]:
        elem_size *= x
    recvcounts = []
    displs = []
    displ = 0
    for r in range(num_node):
        r_start, r_end = get_distributed_range(total_size, r, num_node)
        recvcounts.append((r_end - r_start) * elem_size)
        displs.append(displ)
        displ += recvcounts[-1]
    jk_arr = np.empty((total_size,) + jk_local.shape[1:], dtype=jk_local.dtype)
    comm.Allgatherv(
        jk_local.reshape(-1),
        (jk_arr.reshape(-1), (recvcounts, displs)),
    )
    return jk_arr

@q.timer
def rjackknife(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    rng_state=None,
    n_rand_sample=None,
    jk_blocking_func=None,
    is_normalizing_rand_sample=False,
    is_apply_rand_sample_jk_idx_blocking_shift=True,
    is_use_old_rand_alg=False,
    eps=1,
    is_sync_node=False,
):
    r"""
    Jackknife-bootstrap hybrid resampling.
    Return ``jk_arr``.
    ``len(jk_arr) == 1 + n_rand_sample``
    distribution of ``jk_arr`` should be similar as the distribution of ``avg``.
    ``r_{i,j} ~ N(0, 1)``\n
    ::\n
        if is_normalizing_rand_sample:
            n_j = \sum_i r_{i,j}^2
            r_{i,j} <- \sqrt{n_rand_sample / n_j} r_{i,j}
        data_list_real = [d for d in data_list if d is not None]
        data_arr = np.array(data_list_real, dtype=dtype)
        avg = average(data_arr)
        len(data_list_real) = n
        jk_arr[0] = avg
        jk_arr[i] = avg + \sum_{j=1}^{n} (-eps/\sqrt{n (n - b(i,j))}) r_{i,j} (data_list_real[j] - avg)\n
    where ``b(i,j)`` represent the ``block_size``.\n
    if ``jk_blocking_func`` is provided::\n
        ``jk_blocking_func(i, jk_idx) => blocked jk_idx``\n
    ::\n
        jk_arr[i] = avg + \sum_{j=1}^{n} r_{i,jk_block_func(j)} (jk_arr[j] - avg)\n
    If ``is_sync_node`` is True:\n
        Assume this is a collective operation in a MPI program where every
        node have the same input.  The operation is performed by
        ``rjackknife_sync_node``, which splits the input between the nodes and
        calls ``rjackknife_distributed``; every node obtains the complete
        ``jk_arr``, which agrees with the one obtained with
        ``is_sync_node=False`` up to the floating-point roundoff (not
        bit-for-bit, because the average and the sums over the data set are
        reduced across the nodes).  ``qlat`` (used for ``q.get_comm()``) and
        ``mpi4py`` are imported only when ``is_sync_node`` is True.
    """
    if is_sync_node:
        # Collective MPI operation: every node has the same input and obtains
        # the complete result, which is computed distributedly.
        return rjackknife_sync_node(
            data_list,
            jk_idx_list,
            avg=avg,
            rng_state=rng_state,
            n_rand_sample=n_rand_sample,
            jk_blocking_func=jk_blocking_func,
            is_normalizing_rand_sample=is_normalizing_rand_sample,
            is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
            is_use_old_rand_alg=is_use_old_rand_alg,
            eps=eps,
        )
    if n_rand_sample is None:
        n_rand_sample = 1024
    if rng_state is None:
        rng_state = q.RngState("rejk")
    assert len(data_list) == len(jk_idx_list)
    assert isinstance(n_rand_sample, q.int_types)
    assert n_rand_sample >= 0
    assert isinstance(rng_state, q.RngState)
    if isinstance(data_list, np.ndarray):
        dtype = data_list.dtype
    else:
        dtype = None
    data_list_real = [d for d in data_list if d is not None]
    data_arr = np.array(data_list_real, dtype=dtype)
    if avg is None:
        avg = q.average(data_arr)
    dtype = data_arr.dtype
    jk_idx_list = [jk_idx for jk_idx, d in zip(jk_idx_list, data_list) if d is not None]
    n = len(data_arr)
    #
    r_arr, b_arr = mk_r_i_j_mat(
        n_rand_sample,
        jk_idx_list,
        rng_state,
        jk_blocking_func=jk_blocking_func,
        is_normalizing_rand_sample=is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
        is_use_old_rand_alg=is_use_old_rand_alg,
    )
    n_b_arr = n - b_arr
    n_b_arr[n <= b_arr] = 1
    fac_arr = -eps / np.sqrt(n * n_b_arr)
    fac_arr[n <= b_arr] = 0
    fac_r_arr = fac_arr * r_arr
    pad_shape = (1,) * len(data_arr[0].shape)
    fac_r_arr = fac_r_arr.reshape(fac_r_arr.shape + pad_shape)
    data_diff = data_arr - avg
    jk_rows = avg + np.sum(fac_r_arr * data_diff, axis=1)
    jk_arr = np.empty(
        (
            1 + n_rand_sample,
            *data_arr[0].shape,
        ),
        dtype=dtype,
    )
    jk_arr[0] = avg
    jk_arr[1:] = jk_rows
    return jk_arr

@q.timer
def rjackknife_distributed(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    rng_state=None,
    n_rand_sample=None,
    jk_blocking_func=None,
    is_normalizing_rand_sample=False,
    is_apply_rand_sample_jk_idx_blocking_shift=True,
    is_use_old_rand_alg=False,
    eps=1,
):
    r"""
    ``rjackknife`` with distributed input and distributed output.\n
    This is a collective MPI operation: every node must call it with the same
    parameters and with its own disjoint part of the data set (``data_list``
    and ``jk_idx_list`` are the local parts); every node returns its own part
    of the ``jk_arr`` of ``rjackknife``.  Node ``r`` out of ``num_node`` nodes
    holds the samples in
    ``range(*get_distributed_range(1 + n_rand_sample, r, num_node))``, so
    that::
        jk_arr = np.concatenate(q.get_comm().allgather(jk_local))
    is the complete data set, in the same order as the one obtained with
    ``rjackknife``.  No node needs to hold the whole data set, the whole
    random matrix or the whole result; only the small ``jk_idx`` metadata is
    gathered on every node.\n
    The result agrees with the one of ``rjackknife`` up to the floating-point
    roundoff, but not bit-for-bit: the average and the sums over the data set
    are reduced across the nodes, which changes the order of the
    floating-point additions.\n
    ``rjackknife_sync_node`` performs the corresponding collective operation
    for the case where every node has the whole input; ``g_mk_jk_distributed``
    dispatches to this function.
    """
    fname = "rjackknife_distributed"
    comm, id_node, num_node = get_collective_comm(fname)
    if n_rand_sample is None:
        n_rand_sample = 1024
    if rng_state is None:
        rng_state = q.RngState("rejk")
    (
        data_arr,
        jk_idx_list_local,
        jk_idx_list_glb,
        n,
        elem_shape,
        dtype,
    ) = get_distributed_jk_input(fname, comm, data_list, jk_idx_list)
    avg = get_distributed_avg(comm, data_arr, n, avg)
    #
    total_size = 1 + n_rand_sample
    i_start, i_end = get_distributed_range(total_size, id_node, num_node)
    # ``partial_arr`` holds the contribution of the local data to every sample
    partial_arr = np.zeros((total_size,) + elem_shape, dtype=dtype)
    if i_start == 0 and i_end > i_start:
        # the sample 0 is the average and must be counted only once
        partial_arr[0] = avg
    #
    r_arr, b_arr = mk_r_i_j_mat(
        n_rand_sample,
        jk_idx_list_local,
        rng_state,
        jk_blocking_func=jk_blocking_func,
        is_normalizing_rand_sample=is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
        is_use_old_rand_alg=is_use_old_rand_alg,
        jk_idx_list_for_count=jk_idx_list_glb,
    )
    n_b_arr = n - b_arr
    n_b_arr[n <= b_arr] = 1
    fac_arr = -eps / np.sqrt(n * n_b_arr)
    fac_arr[n <= b_arr] = 0
    fac_r_arr = fac_arr * r_arr
    pad_shape = (1,) * len(elem_shape)
    fac_r_arr = fac_r_arr.reshape(fac_r_arr.shape + pad_shape)
    data_diff = data_arr - avg
    partial_arr[1:] = np.sum(fac_r_arr * data_diff, axis=1)
    return get_reduce_scattered_jk_arr(partial_arr, comm, id_node, num_node, avg)

@q.timer
def rjackknife_sync_node(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    rng_state=None,
    n_rand_sample=None,
    jk_blocking_func=None,
    is_normalizing_rand_sample=False,
    is_apply_rand_sample_jk_idx_blocking_shift=True,
    is_use_old_rand_alg=False,
    eps=1,
):
    r"""
    ``rjackknife`` as a collective MPI operation where every node has the same
    (whole) input and obtains the same complete result.\n
    This implements the ``is_sync_node=True`` option of ``rjackknife``: the
    input is split between the nodes in the order of the nodes and
    ``rjackknife_distributed`` is called on the local parts; the local results
    are then gathered with ``Allgatherv``, so that every node obtains the
    complete ``jk_arr`` of ``rjackknife``.  The result agrees with the one of
    ``rjackknife`` up to the floating-point roundoff, but not bit-for-bit.\n
    ``g_mk_jk_sync_node`` dispatches to this function.
    """
    comm, id_node, num_node = get_collective_comm("rjackknife_sync_node")
    i_start, i_end = get_distributed_range(len(data_list), id_node, num_node)
    jk_local = rjackknife_distributed(
        data_list[i_start:i_end],
        jk_idx_list[i_start:i_end],
        avg=avg,
        rng_state=rng_state,
        n_rand_sample=n_rand_sample,
        jk_blocking_func=jk_blocking_func,
        is_normalizing_rand_sample=is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
        is_use_old_rand_alg=is_use_old_rand_alg,
        eps=eps,
    )
    return get_gathered_jk_arr(jk_local, comm, num_node)

@q.timer
def rjk_mk_jk_val(
    rs_tag,
    val,
    err,
    *,
    n_rand_sample=None,
    rng_state=None,
    eps=1,
):
    """
    return jk_arr
    n = n_rand_sample
    len(jk_arr) == 1 + n
    jk_arr[i] = val + err * r[i] for i in 1..n
    where r[i] ~ N(0, 1)
    """
    if n_rand_sample is None:
        n_rand_sample = 1024
    if rng_state is None:
        rng_state = q.RngState("rejk")
    assert n_rand_sample >= 0
    assert isinstance(rng_state, q.RngState)
    assert isinstance(val, q.real_types)
    assert isinstance(err, q.real_types)
    rs = rng_state.split(str(rs_tag))
    jk_arr = np.zeros((n_rand_sample + 1,), dtype=np.float64)
    jk_arr[0] = val
    r_arr = rs.g_rand_arr((n_rand_sample,))
    r_arr_qnorm = q.qnorm(r_arr)
    r_arr = r_arr * np.sqrt(n_rand_sample / r_arr_qnorm)
    assert abs(q.qnorm(r_arr) / n_rand_sample - 1) < 1e-8
    jk_arr[1:] = val + eps * r_arr * err
    return jk_arr

def rjk_avg(jk_arr):
    return jk_avg(jk_arr)

def rjk_err(jk_arr, eps=1):
    r"""
    Return
    $$
    \frac{1}{eps} \sqrt{ 1/N \sum_{i=1}^N (jk[i] - jk_avg)^2 }.
    $$
    Note: ``
    len(jk_arr) = N + 1.
    jk_avg = jk_arr[0]
    ``
    Same ``eps`` as the ``eps`` used in the ``jackknife`` function.
    """
    avg = jk_avg(jk_arr)
    n = len(jk_arr) - 1
    if n <= 0:
        fac = 1 / abs(eps)
        val = fac * avg
        val = q.filter_np_results(val)
        return val
    diff_sqr = q.average([q.fsqr(jk - avg) for jk in jk_arr[1:]])
    fac = 1 / abs(eps)
    val = fac * q.fsqrt(diff_sqr)
    val = q.filter_np_results(val)
    return val

def rjk_avg_err(rjk_list, eps=1):
    return rjk_avg(rjk_list), rjk_err(rjk_list, eps)

# ----------

default_g_jk_kwargs = dict()

def mk_g_jk_kwargs():
    """
    Return the predefined ``default_g_jk_kwargs``.
    """
    g_jk_kwargs = dict()
    #
    g_jk_kwargs["jk_type"] = "rjk"  # choices: "rjk", "super"
    g_jk_kwargs["eps"] = 1
    #
    # for jk_type = "rjk"
    g_jk_kwargs["n_rand_sample"] = 1024
    g_jk_kwargs["is_normalizing_rand_sample"] = False
    g_jk_kwargs["is_apply_rand_sample_jk_idx_blocking_shift"] = True
    #
    # for jk_type = "super"
    g_jk_kwargs["is_hash_jk_idx"] = True
    g_jk_kwargs["jk_idx_hash_size"] = 1024
    #
    # Is only needed to reproduce old results
    # Possible choice: "v1" (also need default_g_jk_kwargs["is_normalizing_rand_sample"] == False)
    g_jk_kwargs["is_use_old_rand_alg"] = False
    #
    # these parameters are used in jk_blocking_func_default
    g_jk_kwargs["block_size"] = 1
    g_jk_kwargs["block_size_dict"] = {
        "job_tag": 1,
    }
    #
    # Below are items which are not touched in
    # ``get_jk_state`` or ``set_jk_state``
    #
    g_jk_kwargs["rng_state"] = q.RngState("rejk")
    #
    g_jk_kwargs["all_jk_idx"] = None
    g_jk_kwargs["get_all_jk_idx"] = None
    #
    g_jk_kwargs["all_jk_idx_set"] = set()
    #
    # ``is_sync_node`` runs ``g_mk_jk`` as a collective MPI operation.  It only
    # changes how the result is computed (the result agrees up to the
    # floating-point roundoff), so it is deliberately not touched in
    # ``get_jk_state`` or ``set_jk_state`` (and hence not part of the cache key).
    g_jk_kwargs["is_sync_node"] = False
    #
    # jk_blocking_func(i, jk_idx) => blocked_jk_idx
    g_jk_kwargs["jk_blocking_func"] = jk_blocking_func_default
    #
    return g_jk_kwargs

def reset_default_g_jk_kwargs():
    default_g_jk_kwargs.clear()
    default_g_jk_kwargs.update(mk_g_jk_kwargs())

@q.use_kwargs(default_g_jk_kwargs)
def get_jk_state(
    *,
    jk_type,
    eps,
    n_rand_sample,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_hash_jk_idx,
    jk_idx_hash_size,
    is_use_old_rand_alg,
    block_size,
    block_size_dict,
    **_kwargs,
):
    """
    Currently only useful if we set::\n
        q.default_g_jk_kwargs["jk_type"] = "rjk" # this is the default now\n
    and::\n
        q.default_g_jk_kwargs["jk_blocking_func"] = jk_blocking_func_default\n
    Used for ``q.cache_call``.\n
    Example::\n
        @cache_call(get_state=q.get_jk_state)
        def func(...):
            ...
    """
    return (
        jk_type,
        eps,
        n_rand_sample,
        is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift,
        is_hash_jk_idx,
        jk_idx_hash_size,
        is_use_old_rand_alg,
        block_size,
        block_size_dict,
    )

def set_jk_state(state):
    (
        jk_type,
        eps,
        n_rand_sample,
        is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift,
        is_hash_jk_idx,
        jk_idx_hash_size,
        is_use_old_rand_alg,
        block_size,
        block_size_dict,
    ) = state
    g_dict = default_g_jk_kwargs
    g_dict["jk_type"] = jk_type
    g_dict["eps"] = eps
    g_dict["n_rand_sample"] = n_rand_sample
    g_dict["is_normalizing_rand_sample"] = is_normalizing_rand_sample
    g_dict["is_apply_rand_sample_jk_idx_blocking_shift"] = (
        is_apply_rand_sample_jk_idx_blocking_shift
    )
    g_dict["is_hash_jk_idx"] = is_hash_jk_idx
    g_dict["jk_idx_hash_size"] = jk_idx_hash_size
    g_dict["is_use_old_rand_alg"] = is_use_old_rand_alg
    g_dict["block_size"] = block_size
    g_dict["block_size_dict"] = block_size_dict

jk_blocking_traj_shift_arr = q.RngState("jk_blocking_traj_shift_arr").rand_arr(
    16 * 1024
) % (1024 * 1024 * 1024 * 1024)

@q.use_kwargs(default_g_jk_kwargs)
def jk_blocking_func_default(
    i,
    jk_idx,
    *,
    block_size,
    block_size_dict,
    all_jk_idx_set,
    **_kwargs,
):
    """
    return ``blocked_jk_idx``.
    ``blocked_jk_idx`` should uniquely identify the block
    that configuration identified by ``jk_idx`` belongs to.
    The block scheme can be different for different J-B hybrid sample.
    The J-B hybrid sample is indexed by ``i`` (``1 <= i <= n_rand_sample``).
    ``
    block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
    ``
    use default_g_jk_kwargs for
    block_size, block_size_dict, all_jk_idx_set
    """
    if i == 0:
        shift = 0
    else:
        assert i >= 1
        shift = int(
            jk_blocking_traj_shift_arr[(i - 1) % len(jk_blocking_traj_shift_arr)]
        )
    if block_size_dict is None:
        block_size_dict = dict()
    if all_jk_idx_set is not None:
        all_jk_idx_set.add(jk_idx)
    if isinstance(jk_idx, q.int_types):
        traj = jk_idx
        b_shift = shift % block_size
        return (traj + b_shift) // block_size
    elif (
        isinstance(jk_idx, tuple)
        and len(jk_idx) == 2
        and isinstance(jk_idx[1], q.int_types)
    ):
        job_tag, traj = jk_idx
        assert isinstance(job_tag, str)
        assert isinstance(traj, q.int_types)
        block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
        assert isinstance(block_size_for_this_job_tag, q.int_types)
        b_shift = shift % block_size_for_this_job_tag
        return (
            job_tag,
            (traj + b_shift) // block_size_for_this_job_tag,
        )
    else:
        return jk_idx
    assert False

@q.use_kwargs(default_g_jk_kwargs)
@q.timer
def g_mk_jk(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    rng_state,
    jk_blocking_func,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_use_old_rand_alg,
    is_hash_jk_idx,
    jk_idx_hash_size,
    eps,
    is_sync_node=False,
    **_kwargs,
):
    """
    Create a (randomized) Super-Jackknife data set from un-jackknifed data.\n
    ``jk_arr[0]`` is the average of the data and ``jk_arr[1:]`` are the
    resampled values, from which the error is estimated with ``g_jk_avg_err``.
    The data set has ``g_jk_size()`` samples, i.e. ``1 + n_rand_sample`` for
    ``jk_type == "rjk"`` and ``1 + len(all_jk_idx)`` for
    ``jk_type == "super"``, and its dtype is the dtype of the data.\n
    :param data_list: the un-jackknifed data, a list or ``np.ndarray`` of
        values (each value being a ``float``, a ``complex`` or an
        ``np.ndarray``); ``None`` entries are ignored.  For the collective MPI
        variants below the data must have a dtype supported by MPI.
    :param jk_idx_list: the indices that name the entries of ``data_list``,
        with ``len(jk_idx_list) == len(data_list)``, usually
        ``jk_idx_list = [(job_tag, traj,) for traj in traj_list]``.  The
        indices are mapped to the jackknife blocks by ``jk_blocking_func``
        (see ``jk_blocking_func_default``, ``block_size`` and
        ``block_size_dict``).
    :param avg: the average of the data; when ``None`` (the default) it is
        computed from ``data_list``.  Pass a precomputed value to reuse it (it
        must be the average of the whole data set).
    :param is_sync_node: when ``True``, the operation is performed by
        ``g_mk_jk_sync_node``, i.e. as a collective MPI operation in which
        every node holds the whole data set, every node must call this function
        with the same parameters, and every node obtains the complete
        ``jk_arr``.  Use ``g_mk_jk_distributed`` instead when the data set
        itself is split between the nodes.  Both ``jk_type`` values are
        supported; the result agrees with the ``is_sync_node=False`` result up
        to the floating-point roundoff.
    :return: the (randomized) Super-Jackknife data set ``jk_arr``.\n
    The other keyword parameters are the entries of ``default_g_jk_kwargs``,
    which supplies their defaults; set them there, pass them explicitly or use
    the ``q.JkKwargs(...)`` context manager.  The most commonly used are:\n
    - ``jk_type``: ``"rjk"`` (the default) or ``"super"``.
    - ``eps`` (default ``1``): the overall scale of the fluctuations; when the
      data is already jackknifed, multiply it by ``len(data_list)``.
    - ``n_rand_sample`` (default ``1024``): the number of random samples of
      ``"rjk"``.
    - ``is_normalizing_rand_sample``,
      ``is_apply_rand_sample_jk_idx_blocking_shift`` and
      ``is_use_old_rand_alg``: options of the random numbers of ``"rjk"``.
    - ``is_hash_jk_idx``, ``jk_idx_hash_size``, ``all_jk_idx`` and
      ``get_all_jk_idx``: the samples of ``"super"``.
    - ``block_size``, ``block_size_dict`` and ``jk_blocking_func``: the
      jackknife blocks.
    - ``rng_state``: the random numbers of ``"rjk"``.\n
    See ``rjackknife`` and ``sjackknife`` for the formulas, and
    ``docs/qlat-utils/qlat_data.md`` for the full documentation.\n
    Example::
        jk_arr = q.g_mk_jk(data_list, jk_idx_list)
        avg, err = q.g_jk_avg_err(jk_arr)
    """
    if is_sync_node:
        # Collective MPI operation where every node has the same input: only a
        # part of the result is computed on each node and the parts are then
        # gathered, so that every node obtains the complete result.
        jk_arr = g_mk_jk_sync_node(
            data_list,
            jk_idx_list,
            avg=avg,
            jk_type=jk_type,
            all_jk_idx=all_jk_idx,
            get_all_jk_idx=get_all_jk_idx,
            n_rand_sample=n_rand_sample,
            rng_state=rng_state,
            jk_blocking_func=jk_blocking_func,
            is_normalizing_rand_sample=is_normalizing_rand_sample,
            is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
            is_use_old_rand_alg=is_use_old_rand_alg,
            is_hash_jk_idx=is_hash_jk_idx,
            jk_idx_hash_size=jk_idx_hash_size,
            eps=eps,
        )
        return jk_arr
    if jk_type == "super":
        jk_arr = sjackknife(
            data_list,
            jk_idx_list,
            avg=avg,
            is_hash_jk_idx=is_hash_jk_idx,
            jk_idx_hash_size=jk_idx_hash_size,
            rng_state=rng_state,
            all_jk_idx=all_jk_idx,
            get_all_jk_idx=get_all_jk_idx,
            jk_blocking_func=jk_blocking_func,
            eps=eps,
        )
    elif jk_type == "rjk":
        jk_arr = rjackknife(
            data_list,
            jk_idx_list,
            avg=avg,
            n_rand_sample=n_rand_sample,
            rng_state=rng_state,
            jk_blocking_func=jk_blocking_func,
            is_normalizing_rand_sample=is_normalizing_rand_sample,
            is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
            is_use_old_rand_alg=is_use_old_rand_alg,
            eps=eps,
            is_sync_node=is_sync_node,
        )
    else:
        assert False
    return jk_arr

@q.use_kwargs(default_g_jk_kwargs)
@q.timer
def g_mk_jk_distributed(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    rng_state,
    jk_blocking_func,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_use_old_rand_alg,
    is_hash_jk_idx,
    jk_idx_hash_size,
    eps,
    **_kwargs,
):
    """
    Create a (randomized) Super-Jackknife data set when the data set itself is
    split between the MPI nodes.\n
    This is a collective MPI operation: every node must call it with the same
    parameters and with its own disjoint part of the data set, and every node
    returns its own part of the result.  The split of the data set between the
    nodes is free as long as the local parts are disjoint and cover the whole
    data set exactly once (``get_distributed_range(len(data_list), ...)`` is a
    convenient balanced split); no node needs to hold the whole data set, the
    whole random matrix or the whole result, and only the small ``jk_idx``
    metadata is gathered on every node.\n
    The returned ``jk_arr`` is the local part of the data set, which has
    ``g_jk_size()`` samples in total.  The samples are split between the nodes
    in the order of the nodes, independently of the split of the input: node
    ``r`` out of ``num_node`` nodes holds the samples in
    ``range(*get_distributed_range(g_jk_size(...), r, num_node))`` (a node
    holds 0 samples when there are more nodes than samples), so that::
        jk_arr = np.concatenate(q.get_comm().allgather(jk_local))
    is the complete data set, in the same order as the one obtained with
    ``g_mk_jk``.\n
    ``jk_type == "rjk"`` calls ``rjackknife_distributed`` (the random samples
    are split between the nodes and every node computes the contribution of
    its own data to all of them); ``jk_type == "super"`` calls
    ``sjackknife_distributed`` (the samples are the ``all_jk_idx`` entries, or
    the hash based samples).\n
    Requires qlat to be initialized on the whole MPI communicator, i.e. with
    ``q.begin_with_mpi()``, ``q.begin_with_gpt()`` or ``q.begin_with_grid()``
    (or ``q.set_comm(...)``), and the data must have a dtype supported by MPI.\n
    The result agrees with the one of ``g_mk_jk`` (and of
    ``g_mk_jk(..., is_sync_node=True)``) up to the floating-point roundoff, but
    not bit-for-bit: the average and the sums over the data set are reduced
    across the nodes, which changes the order of the floating-point additions.
    The ``is_sync_node`` entry of ``default_g_jk_kwargs`` is ignored, since the
    output of this function is always distributed.\n
    :param data_list: the local part of the un-jackknifed data.
    :param jk_idx_list: the indices that name the local ``data_list``, with
        ``len(jk_idx_list) == len(data_list)``.
    :param avg: the average of the whole data set, the same on every node; when
        ``None`` (the default) it is computed from the data set by summing the
        local contributions over the nodes with ``Allreduce``.
    :return: the local part of the (randomized) Super-Jackknife data set (an
        array with 0 samples on a node which owns no sample).\n
    The other keyword parameters, including ``jk_type`` and ``eps``, are the
    same as for ``g_mk_jk``; their defaults are the entries of
    ``default_g_jk_kwargs``.  See ``g_mk_jk`` for the description of the
    entries.\n
    Example::
        comm = q.get_comm()
        i_start, i_end = q.get_distributed_range(len(data_list), comm.rank, comm.size)
        jk_local = q.g_mk_jk_distributed(
            data_list[i_start:i_end], jk_idx_list[i_start:i_end],
        )
        jk_arr = np.concatenate(comm.allgather(jk_local))
    """
    if jk_type == "super":
        jk_arr = sjackknife_distributed(
            data_list,
            jk_idx_list,
            avg=avg,
            is_hash_jk_idx=is_hash_jk_idx,
            jk_idx_hash_size=jk_idx_hash_size,
            rng_state=rng_state,
            all_jk_idx=all_jk_idx,
            get_all_jk_idx=get_all_jk_idx,
            jk_blocking_func=jk_blocking_func,
            eps=eps,
        )
    elif jk_type == "rjk":
        jk_arr = rjackknife_distributed(
            data_list,
            jk_idx_list,
            avg=avg,
            rng_state=rng_state,
            n_rand_sample=n_rand_sample,
            jk_blocking_func=jk_blocking_func,
            is_normalizing_rand_sample=is_normalizing_rand_sample,
            is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
            is_use_old_rand_alg=is_use_old_rand_alg,
            eps=eps,
        )
    else:
        assert False
    return jk_arr

@q.use_kwargs(default_g_jk_kwargs)
@q.timer
def g_mk_jk_sync_node(
    data_list,
    jk_idx_list,
    *,
    avg=None,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    rng_state,
    jk_blocking_func,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_use_old_rand_alg,
    is_hash_jk_idx,
    jk_idx_hash_size,
    eps,
    **_kwargs,
):
    """
    Create a (randomized) Super-Jackknife data set as a collective MPI
    operation in which every node holds the whole input and obtains the whole
    result.  This is what ``g_mk_jk(..., is_sync_node=True)`` calls.\n
    Every node must call this function with the same parameters and with the
    whole data set (``data_list`` and ``jk_idx_list`` are the complete lists);
    every node then obtains the complete ``jk_arr``, in the same order as the
    one obtained with ``g_mk_jk`` on a single node.  The input is split between
    the nodes in the order of the nodes with ``get_distributed_range`` and the
    local parts are computed by ``rjackknife_sync_node`` or
    ``sjackknife_sync_node`` (which call ``rjackknife_distributed`` or
    ``sjackknife_distributed``) and then gathered, so the work is parallelized
    over the nodes even though every node starts with the whole data set.\n
    Both ``jk_type == "rjk"`` and ``jk_type == "super"`` are supported.  The
    result agrees with the one of ``g_mk_jk`` up to the floating-point
    roundoff, but not bit-for-bit: the average and the sums over the data set
    are reduced across the nodes, which changes the order of the floating-point
    additions.\n
    Requires qlat to be initialized on the whole MPI communicator, i.e. with
    ``q.begin_with_mpi()``, ``q.begin_with_gpt()`` or ``q.begin_with_grid()``
    (or ``q.set_comm(...)``), and the data must have a dtype supported by MPI.\n
    :param data_list: the whole un-jackknifed data, the same on every node.
    :param jk_idx_list: the indices that name the whole ``data_list``, the same
        on every node.
    :param avg: the average of the whole data set, the same on every node; when
        ``None`` (the default) it is computed from the data set by summing the
        local contributions over the nodes with ``Allreduce``.
    :return: the complete (randomized) Super-Jackknife data set ``jk_arr``, the
        same on every node.\n
    The other keyword parameters, including ``jk_type`` and ``eps``, are the
    same as for ``g_mk_jk``; their defaults are the entries of
    ``default_g_jk_kwargs``.  See ``g_mk_jk`` for the description of the
    entries, and ``g_mk_jk_distributed`` for the case in which the data set is
    split between the nodes instead.\n
    Example::
        jk_arr = q.g_mk_jk_sync_node(data_list, jk_idx_list)
        # is equivalent to
        jk_arr = q.g_mk_jk(data_list, jk_idx_list, is_sync_node=True)
    """
    if jk_type == "super":
        jk_arr = sjackknife_sync_node(
            data_list,
            jk_idx_list,
            avg=avg,
            is_hash_jk_idx=is_hash_jk_idx,
            jk_idx_hash_size=jk_idx_hash_size,
            rng_state=rng_state,
            all_jk_idx=all_jk_idx,
            get_all_jk_idx=get_all_jk_idx,
            jk_blocking_func=jk_blocking_func,
            eps=eps,
        )
    elif jk_type == "rjk":
        jk_arr = rjackknife_sync_node(
            data_list,
            jk_idx_list,
            avg=avg,
            rng_state=rng_state,
            n_rand_sample=n_rand_sample,
            jk_blocking_func=jk_blocking_func,
            is_normalizing_rand_sample=is_normalizing_rand_sample,
            is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
            is_use_old_rand_alg=is_use_old_rand_alg,
            eps=eps,
        )
    else:
        assert False
    return jk_arr

@q.use_kwargs(default_g_jk_kwargs)
@q.timer
def g_mk_jk_val(
    rs_tag,
    val,
    err,
    *,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    rng_state,
    is_hash_jk_idx,
    jk_idx_hash_size,
    eps,
    **_kwargs,
):
    """
    Create a jackknife sample with random numbers based on central value ``val`` and error ``err``.\n
    Need::\n
        default_g_jk_kwargs["jk_type"] = "rjk"
        default_g_jk_kwargs["n_rand_sample"] = n_rand_sample
        # e.g. n_rand_sample = 1024
        default_g_jk_kwargs["rng_state"] = rng_state
        # e.g. rng_state = q.RngState("rejk")
    """
    if jk_type == "super":
        jk_val = sjk_mk_jk_val(
            rs_tag,
            val,
            err,
            is_hash_jk_idx=is_hash_jk_idx,
            jk_idx_hash_size=jk_idx_hash_size,
            rng_state=rng_state,
            all_jk_idx=all_jk_idx,
            get_all_jk_idx=get_all_jk_idx,
            eps=eps,
        )
    elif jk_type == "rjk":
        jk_val = rjk_mk_jk_val(
            rs_tag,
            val,
            err,
            n_rand_sample=n_rand_sample,
            rng_state=rng_state,
            eps=eps,
        )
    else:
        assert False
    return jk_val

def g_jk_avg(jk_arr, **_kwargs):
    """
    Return ``avg`` of the ``jk_arr``.
    """
    if isinstance(jk_arr, q.number_types):
        return jk_arr
    return jk_avg(jk_arr)

@q.use_kwargs(default_g_jk_kwargs)
def g_jk_err(jk_arr, *, eps, jk_type, **_kwargs):
    """
    Return ``err`` of the ``jk_arr``.
    """
    if isinstance(jk_arr, q.number_types):
        return 0
    if jk_type == "super":
        return sjk_err(jk_arr, eps=eps)
    elif jk_type == "rjk":
        return rjk_err(jk_arr, eps=eps)
    else:
        assert False
    return None

@q.timer
def g_jk_avg_err(jk_arr, **kwargs):
    """
    Return ``(avg, err,)`` of the ``jk_arr``.
    """
    return g_jk_avg(jk_arr), g_jk_err(jk_arr, **kwargs)

@q.timer
def g_jk_avg_err_arr(jk_arr, **kwargs):
    """
    Return ``avg_err_arr`` of the ``jk_arr``.
    ``
    avg_err_arr.shape = jk_arr[0].shape + (2,)
    ``
    """
    avg, err = g_jk_avg_err(jk_arr, **kwargs)
    avg_err_arr = np.stack(
        [
            avg,
            err,
        ]
    )
    avg_err_arr = np.moveaxis(avg_err_arr, 0, -1).copy()
    return avg_err_arr

@q.use_kwargs(default_g_jk_kwargs)
def g_jk_size(
    *,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    is_hash_jk_idx,
    jk_idx_hash_size,
    **_kwargs,
):
    """
    Return number of samples for the (randomized) Super-Jackknife data set.
    """
    if jk_type == "super":
        if all_jk_idx is None:
            if get_all_jk_idx is None:
                assert is_hash_jk_idx
                all_jk_idx = [
                    "avg",
                ] + list(range(jk_idx_hash_size))
            else:
                all_jk_idx = get_all_jk_idx()
        assert all_jk_idx[0] == "avg"
        n_super_sample = len(all_jk_idx) - 1
        assert n_super_sample >= 0
        return 1 + n_super_sample
    elif jk_type == "rjk":
        return 1 + n_rand_sample
    else:
        assert False
    return None

@q.use_kwargs(default_g_jk_kwargs)
def g_jk_blocking_func(
    i,
    jk_idx,
    *,
    jk_blocking_func,
    **_kwargs,
):
    """
    Return ``jk_blocking_func(jk_idx)``.
    """
    if jk_blocking_func is None:
        return jk_idx
    else:
        return jk_blocking_func(i, jk_idx)

@q.use_kwargs(default_g_jk_kwargs)
def g_jk_sample_size(
    job_tag,
    traj_list,
    **_kwargs,
):
    jk_idx_list = [
        (
            job_tag,
            traj,
        )
        for traj in traj_list
    ]
    b_jk_idx_set = set(
        g_jk_blocking_func(0, jk_idx, **kwargs) for jk_idx in jk_idx_list
    )
    return len(b_jk_idx_set)

reset_default_g_jk_kwargs()

# ----

class JkKwargs(q.NewDictValues):
    """
    Example:
    #
    with q.JkKwargs(n_rand_sample=1024, block_size=10, block_size_dict={ "48I": 20, }):
        ...
    #
    """

    def __init__(self, **kwargs):
        super().__init__(default_g_jk_kwargs, **kwargs)

# ----

# ---- old funcs

def merge_jk_idx(*jk_idx_list_list):
    for jk_idx_list in jk_idx_list_list:
        assert jk_idx_list[0] == "avg"
    return [
        "avg",
    ] + [jk_idx for jk_idx_list in jk_idx_list_list for jk_idx in jk_idx_list[1:]]

@q.timer
def rejk_list(jk_list, jk_idx_list, all_jk_idx):
    """
    Super jackknife
    ``jk_idx_list`` should be contained in ``all_jk_idx`` and have the same order.
    Does not properly honor the (N-1) formula in error calculation.
    """
    assert jk_idx_list[0] == "avg"
    assert all_jk_idx[0] == "avg"
    assert len(jk_idx_list) == len(jk_list)
    assert len(jk_idx_list) <= len(all_jk_idx)
    is_np_arr = isinstance(jk_list, np.ndarray)
    jk_avg = jk_list[0]
    size_new = len(all_jk_idx)
    i_new = 0
    jk_list_new = []
    for i, idx in enumerate(jk_idx_list):
        while all_jk_idx[i_new] != idx:
            jk_list_new.append(jk_avg)
            i_new += 1
            assert i_new < size_new
        jk_list_new.append(jk_list[i])
        i_new += 1
    while i_new < size_new:
        jk_list_new.append(jk_avg)
        i_new += 1
    assert i_new == size_new
    assert size_new == len(jk_list_new)
    if is_np_arr:
        jk_list_new = np.array(jk_list_new, dtype=jk_list.dtype)
    return jk_list_new

@q.timer
def rjk_jk_list(
    jk_list,
    jk_idx_list,
    n_rand_sample,
    rng_state,
    jk_blocking_func=None,
    is_normalizing_rand_sample=False,
    is_apply_rand_sample_jk_idx_blocking_shift=True,
    is_use_old_rand_alg=False,
):
    r"""
    return jk_list
    len(jk_list) == 1 + n_rand_sample
    distribution of jk_list should be similar as the distribution of avg
    r_{i,j} ~ N(0, 1)
    if is_normalizing_rand_sample:
        n_j = \sum_i r_{i,j}^2
        r_{i,j} <- \sqrt{n_rand_sample / n_j} r_{i,j}
    avg = jk_list[0]
    len(jk_list) = n + 1
    jk_list[i] = avg + \sum_{j=1}^{n} r_{i,j} (jk_list[j] - avg)
    #
    if ``jk_blocking_func`` is provided:
    ``
    jk_blocking_func(i, jk_idx) => blocked jk_idx
    ``
    Note that: ``1 <= i <= n_rand_sample``
    ``
    jk_list[i] = avg + \sum_{j=1}^{n} r_{i,jk_block_func(i, j)} (jk_list[j] - avg)
    ``
    """
    assert jk_idx_list[0] == "avg"
    assert isinstance(n_rand_sample, q.int_types)
    assert n_rand_sample >= 0
    assert isinstance(rng_state, q.RngState)
    is_np_arr = isinstance(jk_list, np.ndarray)
    n = len(jk_list) - 1
    r_arr, b_arr = mk_r_i_j_mat(
        n_rand_sample,
        jk_idx_list[1:],
        rng_state,
        jk_blocking_func=jk_blocking_func,
        is_normalizing_rand_sample=is_normalizing_rand_sample,
        is_apply_rand_sample_jk_idx_blocking_shift=is_apply_rand_sample_jk_idx_blocking_shift,
        is_use_old_rand_alg=is_use_old_rand_alg,
    )
    avg = jk_list[0]
    if is_np_arr:
        jk_arr = jk_list
        jk_diff = jk_arr[1:] - avg
        rjk_arr = np.empty(
            (
                1 + n_rand_sample,
                *avg.shape,
            ),
            dtype=jk_arr.dtype,
        )
        rjk_arr[:] = avg
        for j in range(n):
            for i in range(n_rand_sample):
                rjk_arr[i + 1] += r_arr[i, j] * jk_diff[j]
        return rjk_arr
    else:
        rjk_list = [
            avg,
        ]
        jk_diff = [jk_list[j] - avg for j in range(1, n + 1)]
        for i in range(n_rand_sample):
            rjk_list.append(avg + sum([r_arr[i, j] * jk_diff[j] for j in range(n)]))
        return rjk_list

@q.use_kwargs(default_g_jk_kwargs)
@q.timer
def g_jk(data_list, *, eps, **_kwargs):
    """
    Obsolete, call ``g_mk_jk`` instead.
        --
    Perform initial Jackknife for the original data set.\n
    """
    return jackknife(data_list, eps=eps)

@q.use_kwargs(default_g_jk_kwargs)
@q.timer
def g_rejk(
    jk_list,
    jk_idx_list,
    *,
    jk_type,
    all_jk_idx,
    get_all_jk_idx,
    n_rand_sample,
    rng_state,
    jk_blocking_func,
    is_normalizing_rand_sample,
    is_apply_rand_sample_jk_idx_blocking_shift,
    is_use_old_rand_alg,
    **_kwargs,
):
    """
    Obsolete, call ``g_mk_jk`` instead.
        --
    Perform (randomized) Super-Jackknife for the Jackknife data set.
        --
    :jk_list: usually the Jackknife data set obtained with ``g_jk(data_list)``.
    :jk_idx_list: should be list of indices that names the ``jk_list``.
    :jk_type: ``[ "rjk", "super", ]``
    :returns: (randomized) Super-Jackknife data set.
    Note that::
        len(jk_list) == len(jk_idx_list)
        jk_idx_list[0] == "avg"
    """
    if jk_type == "super":
        if jk_blocking_func is not None:
            q.displayln_info(
                f"g_rejk: jk_type={jk_type} does not support jk_blocking_func={jk_blocking_func}"
            )
        if all_jk_idx is None:
            assert get_all_jk_idx is not None
            all_jk_idx = get_all_jk_idx()
        return rejk_list(
            jk_list,
            jk_idx_list,
            all_jk_idx,
        )
    elif jk_type == "rjk":
        return rjk_jk_list(
            jk_list,
            jk_idx_list,
            n_rand_sample,
            rng_state,
            jk_blocking_func,
            is_normalizing_rand_sample,
            is_apply_rand_sample_jk_idx_blocking_shift,
            is_use_old_rand_alg,
        )
    else:
        assert False
    return None

# ----

def mk_jk_blocking_func(block_size=1, block_size_dict=None, all_jk_idx_set=None):
    """
    Recommend to use ``jk_blocking_func_default`` instead.
    #
    block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
    """
    if block_size_dict is None:
        block_size_dict = dict()
    #
    def jk_blocking_func(jk_idx):
        if all_jk_idx_set is not None:
            all_jk_idx_set.add(jk_idx)
        if isinstance(jk_idx, q.int_types):
            traj = jk_idx
            return traj // block_size
        elif (
            isinstance(jk_idx, tuple)
            and len(jk_idx) == 2
            and isinstance(jk_idx[1], q.int_types)
        ):
            job_tag, traj = jk_idx
            assert isinstance(job_tag, str)
            assert isinstance(traj, q.int_types)
            block_size_for_this_job_tag = block_size_dict.get(job_tag, block_size)
            assert isinstance(block_size_for_this_job_tag, q.int_types)
            return (
                job_tag,
                traj // block_size_for_this_job_tag,
            )
        else:
            return jk_idx
    #
    return jk_blocking_func

# ----

def add_jk_idx(arr):
    """
    arr: no jk index
    return: add trivial jk index in the LAST axis
    """
    return arr.reshape(arr.shape + (1,))

def jk_transpose(arr):
    """
    arr: jk index is the 0th axis
    return: jk index is the last axis
    """
    shape = arr.shape
    ndim = len(shape)
    if ndim <= 1:
        return arr
    axes = list(range(1, ndim)) + [
        0,
    ]
    return arr.transpose(axes)

def jk_transpose_back(arr):
    """
    jk_transpose_back(jk_transpose(arr)) == arr
    """
    shape = arr.shape
    ndim = len(shape)
    if ndim <= 1:
        return arr
    axes = [
        ndim - 1,
    ] + list(range(0, ndim - 1))
    return arr.transpose(axes)
