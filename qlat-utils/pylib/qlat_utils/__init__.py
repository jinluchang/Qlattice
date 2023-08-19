"""
Qlattice utility package\n
Usage::\n
    import qlat_utils as q\n
Will also be loaded by ``import qlat as q`` together with other ``qlat`` functions.
"""

from . import c

from .c import \
        Coordinate, \
        mod, \
        smod, \
        middle_mod, \
        coordinate_from_index, \
        index_from_coordinate, \
        RngState, \
        as_wilson_matrix, \
        as_wilson_matrix_g5_herm, \
        benchmark_matrix_functions

from .c import \
        WilsonMatrix, \
        SpinMatrix, \
        ColorMatrix, \
        get_gamma_matrix, \
        wilson_matrix_g5_herm, \
        mat_tr_sm, \
        mat_tr_cm, \
        mat_tr_wm, \
        mat_tr_wm_wm, \
        mat_tr_wm_sm, \
        mat_tr_sm_wm, \
        mat_tr_sm_sm, \
        mat_tr_wm_cm, \
        mat_tr_cm_wm, \
        mat_tr_cm_cm, \
        mat_mul_wm_wm, \
        mat_mul_wm_sm, \
        mat_mul_sm_wm, \
        mat_mul_sm_sm, \
        mat_mul_wm_cm, \
        mat_mul_cm_wm, \
        mat_mul_cm_cm

from .c import \
        ElemType, \
        ElemTypeColorMatrix, \
        ElemTypeWilsonMatrix, \
        ElemTypeNonRelWilsonMatrix, \
        ElemTypeIsospinMatrix, \
        ElemTypeSpinMatrix, \
        ElemTypeWilsonVector, \
        ElemTypeComplex, \
        ElemTypeComplexF, \
        ElemTypeDouble, \
        ElemTypeFloat, \
        ElemTypeLong, \
        ElemTypeInt64t, \
        ElemTypeInt8t, \
        ElemTypeChar

from .timer import *

from .ama import *

from .load_prop import *

from .cache import *

from .qar import *

from .utils import *

from .utils_io import *

from .lat_io import *

from .data import *

from .qplot import *

from .parallel import *

from .get_include_dir import *

verbose_level("default")
