from .c import \
        as_wilson_matrix, \
        as_wilson_matrix_g5_herm

from .ama import \
        ama_apply1

def load_prop(x):
    if isinstance(x, tuple):
        assert len(x) == 2 and x[0] == "g5_herm"
        return ama_apply1(as_wilson_matrix_g5_herm, x[1])
    return ama_apply1(as_wilson_matrix, x)
