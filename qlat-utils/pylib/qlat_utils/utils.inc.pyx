def random_permute(list l, RngState rs):
    # Do not change ``l''.
    # Return a new permuted list.
    cdef long size = len(l)
    cdef cc.std_vector[cc.PyObject*] vec
    vec.resize(size)
    cdef long i = 0
    for i in range(size):
        vec[i] = <cc.PyObject*>l[i]
    cc.random_permute(vec, rs.xx)
    cdef list l_new = []
    for i in range(size):
        l_new.append(<object>vec[i])
    return l_new

def displayln_malloc_stats():
    cc.displayln_malloc_stats()

### -------------------------------------------------------------------

def get_all_caches_info():
    cdef list l = cc.get_all_caches_info()
    return [ str(v) for v in l ]

def clear_all_caches():
    cc.clear_all_caches()
