import numpy as np

use_reference_implementation = False

mpi_comm = None

def set_mpi_comm(comm):
    global mpi_comm
    mpi_comm = comm

def get_mpi_comm():
    global mpi_comm
    if mpi_comm is None:
        from mpi4py import MPI
        mpi_comm = MPI.COMM_WORLD
    return mpi_comm

def bcast_py(x, root=0, comm=None):
    if comm is None:
        comm = get_mpi_comm()
    return comm.bcast(x, root)

class DistArray:

    """
    self.n : total size for the first dimension (the distributed dimension) of the array.
    self.x : np.ndarray for this MPI process (pad zeros if does not divide).
    self.comm : comm used when generate this array
    """

    def __init__(self, *, comm=None):
        if comm is None:
            comm = get_mpi_comm()
        self.x = np.zeros(1, dtype=np.float64)
        self.n = 1
        self.comm = comm

    def __repr__(self):
        return f"DistArray({self.n!r},\n{self.x!r})"

    def sum_ref(self, axis=None, *, keepdims=False):
        """
        return the results of sum, never as DistArray.
        Collective operation.
        """
        v = all_gather_arr(self)
        return v.sum(axis=axis, keepdims=keepdims)

    def sum(self, axis=None, *, keepdims=False):
        """
        return the results of sum, never as DistArray.
        Collective operation.
        """
        if use_reference_implementation:
            return self.sum_ref(axis=axis, keepdims=keepdims)
        comm = self.comm
        if axis is not None:
            if isinstance(axis, int):
                axis = (axis,)
            if 0 not in axis:
                d_v = DistArray(comm=comm)
                d_v.comm = comm
                d_v.n = self.n
                d_v.x = self.x.sum(axis=axis, keepdims=keepdims)
                return all_gather_arr(d_v)
        v = self.x.sum(axis=axis, keepdims=True)
        comm.Allreduce(v.copy(), v)
        return v.sum(axis=axis, keepdims=keepdims)

    def __add__(self, v):
        comm = self.comm
        d_ret = DistArray(comm=comm)
        d_ret.n = self.n
        if isinstance(v, DistArray):
            assert self.comm is v.comm
            assert self.n == v.n
            assert len(self.x.shape) == len(v.x.shape)
            d_ret.x = self.x + v.x
        else:
            if isinstance(v, np.ndarray):
                assert len(self.x.shape) >= len(v.shape)
            d_ret.x = self.x + v
        return d_ret

    def __sub__(self, v):
        comm = self.comm
        d_ret = DistArray(comm=comm)
        d_ret.n = self.n
        if isinstance(v, DistArray):
            assert self.comm is v.comm
            assert self.n == v.n
            assert len(self.x.shape) == len(v.x.shape)
            d_ret.x = self.x - v.x
        else:
            if isinstance(v, np.ndarray):
                assert len(self.x.shape) >= len(v.shape)
            d_ret.x = self.x - v
        return d_ret

    def __mul__(self, v):
        comm = self.comm
        d_ret = DistArray(comm=comm)
        d_ret.n = self.n
        if isinstance(v, DistArray):
            assert self.comm is v.comm
            assert self.n == v.n
            assert len(self.x.shape) == len(v.x.shape)
            d_ret.x = self.x * v.x
        else:
            if isinstance(v, np.ndarray):
                assert len(self.x.shape) >= len(v.shape)
            d_ret.x = self.x * v
        return d_ret

    def __truediv__(self, v):
        comm = self.comm
        d_ret = DistArray(comm=comm)
        d_ret.n = self.n
        if isinstance(v, DistArray):
            assert self.comm is v.comm
            assert self.n == v.n
            assert len(self.x.shape) == len(v.x.shape)
            d_ret.x = self.x / v.x
        else:
            if isinstance(v, np.ndarray):
                assert len(self.x.shape) >= len(v.shape)
            d_ret.x = self.x / v
        return d_ret

    __radd__ = __add__

    __rmul__ = __mul__

    def __rsub__(self, v):
        comm = self.comm
        d_ret = DistArray(comm=comm)
        d_ret.n = self.n
        if isinstance(v, DistArray):
            assert self.comm is v.comm
            assert self.n == v.n
            assert len(self.x.shape) == len(v.x.shape)
            d_ret.x = v.x - self.x
        else:
            if isinstance(v, np.ndarray):
                assert len(self.x.shape) >= len(v.shape)
            d_ret.x = v - self.x
        return d_ret

    def __rtruediv__(self, v):
        comm = self.comm
        d_ret = DistArray(comm=comm)
        d_ret.n = self.n
        if isinstance(v, DistArray):
            assert self.comm is v.comm
            assert self.n == v.n
            assert len(self.x.shape) == len(v.x.shape)
            d_ret.x = v.x / self.x
        else:
            if isinstance(v, np.ndarray):
                assert len(self.x.shape) >= len(v.shape)
            d_ret.x = v / self.x
        return d_ret

    def __matmul__(self, v):
        return d_matmul(self, v)

    def transpose2d_ref(self):
        """
        return DistArray
        transpose the first two dims
        Collective operation.
        """
        assert len(self.x.shape) >= 2
        comm = self.comm
        rank = comm.Get_rank()
        root = 0
        axes = list(range(len(self.x.shape)))
        axes[:2] = [ 1, 0, ]
        v = gather_arr(self, root)
        if rank == root:
            v = v.transpose(axes)
        d_v = scatter_arr(v, root, self.comm)
        return d_v

    def transpose2d(self):
        """
        return DistArray
        transpose the first two dims
        Collective operation.
        """
        if use_reference_implementation:
            return self.transpose2d_ref()
        comm = self.comm
        size = comm.Get_size()
        axes = list(range(len(self.x.shape)))
        axes[:2] = [ 1, 0, ]
        vec = self.x.transpose(axes)
        vec_len = len(vec)
        d_vec_len = (vec_len - 1) // size + 1
        new_vec_len = d_vec_len * size
        assert new_vec_len >= vec_len
        if new_vec_len > vec_len:
            # pad zeros to vec
            new_vec = np.empty((new_vec_len,) + vec.shape[1:], dtype=vec.dtype)
            new_vec[:vec_len] = vec
            new_vec[vec_len:] = 0
            vec = new_vec
        else:
            # make sure it is contiguous
            assert d_vec_len * size == vec_len
        vec = np.ascontiguousarray(vec)
        nvec = np.empty((size, d_vec_len,) + vec.shape[1:], dtype=vec.dtype)
        comm.Alltoall(vec, nvec)
        axes = list(range(len(nvec.shape)))
        axes[:2] = [ 1, 0, ]
        nvec = nvec.transpose(axes)
        nvec = nvec.reshape(d_vec_len, size * vec.shape[1], *vec.shape[2:])
        assert nvec.shape[1] >= self.n
        nvec = nvec[:, :self.n]
        nvec = np.ascontiguousarray(nvec)
        d_ret = DistArray(comm=comm)
        d_ret.n = vec_len
        d_ret.x = nvec
        return d_ret

    def transpose(self, axes=None):
        shape = self.x.shape
        assert len(shape) >= 1
        if len(shape) == 1:
            return self
        elif len(shape) == 2:
            if (axes is None) or (tuple(axes) == (1, 0,)):
                return self.transpose2d()
            else:
                assert tuple(axes) == (0, 1,)
                return self
        else:
            raise Exception(f"DistArray.transpose: axes={axes}")

    def conj(self):
        """
        return DistArray
        Collective operation.
        """
        comm = self.comm
        d_ret = DistArray(comm=comm)
        d_ret.n = self.n
        d_ret.x = self.x.conj()
        return d_ret

###

def d_matmul_ref(d_mat, d_vec):
    """
    return DistArray
    d_mat is distributed mat
    d_vec is distributed vec
    """
    assert d_mat.comm is d_vec.comm
    assert len(d_mat.x.shape) >= 2
    comm = d_mat.comm
    rank = comm.Get_rank()
    root = 0
    mat = gather_arr(d_mat, root)
    vec = gather_arr(d_vec, root)
    if root == rank:
        ret = np.matmul(mat, vec)
    else:
        ret = None
    d_ret = scatter_arr(ret, root, comm)
    return d_ret

def d_matmul(d_mat, d_vec):
    """
    return DistArray
    d_mat is distributed mat
    d_vec is distributed vec
    """
    if use_reference_implementation:
        return d_matmul_ref(d_mat, d_vec)
    assert d_mat.comm is d_vec.comm
    assert len(d_mat.x.shape) >= 2
    comm = d_mat.comm
    d_ret = DistArray(comm=comm)
    vec = all_gather_arr(d_vec)
    ret = np.matmul(d_mat.x, vec)
    ret = np.ascontiguousarray(ret)
    d_ret.n = d_mat.n
    d_ret.x = ret
    return d_ret

def d_trace_ref(d_mat):
    assert len(d_mat.x.shape) >= 2
    mat = all_gather_arr(d_mat)
    return np.trace(mat)

def d_trace(d_mat):
    if use_reference_implementation:
        return d_trace_ref(d_mat)
    assert len(d_mat.x.shape) >= 2
    comm = d_mat.comm
    size = comm.Get_size()
    rank = comm.Get_rank()
    d_vec_len = len(d_mat.x)
    r = np.trace(d_mat.x, offset=rank * d_vec_len)
    if isinstance(r, np.ndarray):
        comm.Allreduce(r.copy(), r)
    else:
        r = comm.allreduce(r)
    return r

###

def scatter_arr(vec, root=0, comm=None):
    """
    return DistArray
    only need value of vec on the root node
    Pad zeros to vec if len(vec) % comm.Get_size() != 0
    """
    if comm is None:
        comm = get_mpi_comm()
    d_vec = DistArray(comm=comm)
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank == root:
        assert isinstance(vec, np.ndarray)
        shape_dtype = (vec.shape, vec.dtype,)
    else:
        shape_dtype = None
    shape, dtype = comm.bcast(shape_dtype, root)
    assert len(shape) >= 1
    vec_len = shape[0]
    d_vec_len = (vec_len - 1) // size + 1
    new_vec_len = d_vec_len * size
    assert new_vec_len >= vec_len
    if rank == root:
        if new_vec_len > vec_len:
            # pad zeros to vec
            new_vec = np.empty((new_vec_len,) + shape[1:], dtype=dtype)
            new_vec[:vec_len] = vec
            new_vec[vec_len:] = 0
            vec = new_vec
        # make sure vec is contiguous
        vec = np.ascontiguousarray(vec)
    d_vec.n = vec_len
    d_vec.x = np.empty((d_vec_len,) + shape[1:], dtype=dtype)
    comm.Scatter(vec, d_vec.x, root)
    d_vec.comm = comm
    return d_vec

def gather_arr(d_vec, root=0):
    """
    return np.ndarray on node root
    return None on other nodes
    collective operation, need to run on all nodes.
    Remove padded zeros added
    """
    assert isinstance(d_vec, DistArray)
    comm = d_vec.comm
    size = comm.Get_size()
    rank = comm.Get_rank()
    d_vec_len = len(d_vec.x)
    new_vec_len = d_vec_len * size
    if rank == root:
        vec_g = np.zeros((new_vec_len,) + d_vec.x.shape[1:], dtype=d_vec.x.dtype)
    else:
        vec_g = None
    d_vec.x = np.ascontiguousarray(d_vec.x)
    comm.Gather(d_vec.x, vec_g, root)
    if rank == root:
        assert new_vec_len >= d_vec.n
        vec_g = vec_g[:d_vec.n]
        # vec_g = np.ascontiguousarray(vec_g)
    return vec_g

def all_gather_arr(d_vec):
    """
    return np.ndarray on all nodes
    Remove padded zeros
    """
    assert isinstance(d_vec, DistArray)
    comm = d_vec.comm
    size = comm.Get_size()
    d_vec_len = len(d_vec.x)
    new_vec_len = d_vec_len * size
    vec_g = np.zeros((new_vec_len,) + d_vec.x.shape[1:], dtype=d_vec.x.dtype)
    d_vec.x = np.ascontiguousarray(d_vec.x)
    comm.Allgather(d_vec.x, vec_g)
    assert new_vec_len >= d_vec.n
    vec_g = vec_g[:d_vec.n]
    # vec_g = np.ascontiguousarray(vec_g)
    return vec_g
