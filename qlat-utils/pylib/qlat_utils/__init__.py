"""
Qlattice utility package\n
Usage::\n
    import qlat_utils as q\n
Will also be loaded by ``import qlat as q`` together with other ``qlat`` functions.
"""

import qlat_utils.c as c

from qlat_utils.c import \
        Coordinate, \
        RngState, \
        as_wilson_matrix, \
        as_wilson_matrix_g5_herm, \
        benchmark_matrix_functions

from qlat_utils.c import \
        WilsonMatrix, \
        SpinMatrix, \
        get_gamma_matrix, \
        mat_tr_sm, \
        mat_tr_wm, \
        mat_tr_wm_wm, \
        mat_tr_wm_sm, \
        mat_tr_sm_wm, \
        mat_tr_sm_sm, \
        mat_mul_wm_wm, \
        mat_mul_wm_sm, \
        mat_mul_sm_wm, \
        mat_mul_sm_sm

from qlat_utils.c import \
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

from qlat_utils.timer import *

from qlat_utils.cache import *

from qlat_utils.qar import *

from qlat_utils.utils import *

from qlat_utils.utils_io import *

from qlat_utils.lat_io import *

from qlat_utils.data import *

from qlat_utils.qplot import *

from qlat_utils.parallel import *

from qlat_utils.get_include_dir import *

verbose_level("default")

del c
