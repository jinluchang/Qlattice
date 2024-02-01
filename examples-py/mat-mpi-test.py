#!/usr/bin/env python3

import numpy as np

from qlat.mat_mpi import get_mpi_comm, set_mpi_comm, DistArray, scatter_arr, all_gather_arr, gather_arr, d_matmul, d_trace

comm = get_mpi_comm()
size = comm.Get_size()
rank = comm.Get_rank()

print(f"size={size} ; rank={rank}")

root = 0

if rank == root:
    vec = np.arange(25, dtype=np.complex128) + 1j
    vec = vec.reshape(5, 5)
else:
    vec = None

d_vec = scatter_arr(vec, root)

print(rank, d_vec)
print(rank, all_gather_arr(d_vec))

print(rank, gather_arr(d_vec, root))

d_vec2 = d_matmul(d_vec, d_vec)

vec2 = gather_arr(d_vec2, root)
if rank == root:
    print(vec2)

d_vec2 = d_vec @ d_vec
vec2 = gather_arr(d_vec2, root)
if rank == root:
    print(vec2)

vec2 = all_gather_arr(d_vec) @ all_gather_arr(d_vec)
if rank == root:
    print(vec2)

axis=()
keepdims=False
r = d_vec2.sum(axis, keepdims=keepdims)
if rank == root:
    print()
    print(r)
    print(vec2.sum(axis, keepdims=keepdims))

d_vec2t = d_vec2.transpose2d()
vec2t = gather_arr(d_vec2t, root)
if rank == root:
    print(vec2t)

d_vec3 = 3 * d_vec2t * 2
vec3 = gather_arr(d_vec3, root)
if rank == root:
    print(vec3)

tr = d_trace(d_vec3)
if rank == root:
    print(tr)

print("simple demo")

arr_list = []
nt = 4
for i in range(nt):
    arr = None
    if rank == (i % size):
        arr = (np.arange(16.0) + 1j).reshape(4, 4) # loading data from this tslsice locally
    arr_list.append(arr)

lam = np.arange(1.0, 5.0)
print(lam)

# d_arr[i, j] / lam[j] -> d_arr[i, j]

d_arr_list = []
d_arr_list_div = []
for i in range(nt):
    arr = arr_list[i]
    d_arr = scatter_arr(arr, root=i % size)
    d_arr_list.append(d_arr)
    d_arr_list_div.append(d_arr / lam)

d_arr_list_tt = []

for i in range(nt):
    d_arr_list_tt.append(d_arr_list[i].transpose() / lam)

arr_list = None

# print(rank, d_arr_list[1])

print(rank, (d_trace(d_arr_list[1] @ d_arr_list[0])))

s = 0

for i in range(nt):
    for j in range(nt):
        v = (d_arr_list_div[i] * d_arr_list_tt[j]).sum()
        s += v
print(s / nt**2)

comm.barrier()

if rank == 0:
    print(f"CHECK: finished successfully.")
