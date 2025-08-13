# ----------------------------------

__all__ = []

__all__ += [
        'CoordinateD',
        'Coordinate',
        'mod_coordinate',
        'smod_coordinate',
        'smod_sym_coordinate',
        'middle_mod_coordinate',
        'coordinate_from_index',
        'index_from_coordinate',
        'mod_coordinate_d',
        'smod_coordinate_d',
        'smod_sym_coordinate_d',
        'middle_mod_coordinate_d',
        ]

__all__ += [
        'RngState',
        'get_data_sig',
        'random_permute',
        ]

__all__ += [
        'as_wilson_matrix',
        'as_wilson_matrix_g5_herm',
        'benchmark_matrix_functions',
        'WilsonMatrix',
        'SpinMatrix',
        'ColorMatrix',
        'get_gamma_matrix',
        'wilson_matrix_g5_herm',
        'mat_tr_sm',
        'mat_tr_cm',
        'mat_tr_wm',
        'mat_tr_wm_wm',
        'mat_tr_wm_sm',
        'mat_tr_sm_wm',
        'mat_tr_sm_sm',
        'mat_tr_wm_cm',
        'mat_tr_cm_wm',
        'mat_tr_cm_cm',
        'mat_mul_a_wm',
        'mat_mul_a_sm',
        'mat_mul_a_cm',
        'mat_mul_wm_wm',
        'mat_mul_wm_sm',
        'mat_mul_sm_wm',
        'mat_mul_sm_sm',
        'mat_mul_wm_cm',
        'mat_mul_cm_wm',
        'mat_mul_cm_cm',
        'mat_add_wm_wm',
        'mat_add_sm_sm',
        'mat_add_cm_cm',
        'mat_epsilon_contraction_wm_wm_wm',
        ]

__all__ += [
        'ElemType',
        'ElemTypeColorMatrix',
        'ElemTypeWilsonMatrix',
        'ElemTypeNonRelWilsonMatrix',
        'ElemTypeIsospinMatrix',
        'ElemTypeSpinMatrix',
        'ElemTypeWilsonVector',
        'ElemTypeComplexD',
        'ElemTypeComplexF',
        'ElemTypeRealD',
        'ElemTypeRealF',
        'ElemTypeLong',
        'ElemTypeInt',
        'ElemTypeInt64t',
        'ElemTypeInt32t',
        'ElemTypeInt8t',
        'ElemTypeChar',
        ]

__all__ += [
        'get_id_node',
        'get_num_node',
        'sync_node',
        'get_verbose_level',
        'set_verbose_level',
        'get_time',
        'get_start_time',
        'set_start_time',
        'get_actual_start_time',
        'set_actual_start_time',
        'get_total_time',
        'get_actual_total_time',
        'get_time_limit',
        'set_time_limit',
        'get_remaining_time',
        'get_time_budget',
        'set_time_budget',
        'timer_display',
        'timer_autodisplay',
        'timer_display_stack',
        'timer_display_stack_always',
        'timer_reset',
        'timer_fork',
        'timer_merge',
        'timer',
        'timer_fname',
        'timer_verbose',
        'timer_verbose_fname',
        'timer_flops',
        'timer_flops_fname',
        'timer_verbose_flops',
        'timer_verbose_flops_fname',
        'Timer',
        'TimerNone',
        'displayln',
        'displayln_info',
        'set_display_method',
        ]

for ext, ext2 in [ ('', '',), ('RealF', '_real_f',), ('Int', '_int',), ('Long', '_long',), ]:
    __all__ += [
            f'LatData{ext}',
            f'mk_lat_data{ext2}',
            f'load_lat_data{ext2}',
            ]

__all__ += [
        'basename',
        'dirname',
        'all_dirname_vec',
        'remove_trailing_slashes',
        'qls',
        'qls_all',
        'does_file_exist',
        'is_directory',
        'is_regular_file',
        'qmkdir',
        'qmkdir_p',
        'qrename',
        'qremove',
        'qremove_all',
        'qtruncate',
        #
        'clear_is_directory_cache',
        'remove_entry_directory_cache',
        'is_directory_cache',
        'is_regular_file_cache',
        'does_file_exist_cache',
        #
        'qmkdir_info',
        'qmkdir_p_info',
        'qrename_info',
        'qremove_info',
        'qremove_all_info',
        #
        'qls_sync_node',
        'qls_all_sync_node',
        'does_file_exist_sync_node',
        'is_directory_sync_node',
        'is_regular_file_sync_node',
        'qmkdir_sync_node',
        'qmkdir_p_sync_node',
        'qremove_sync_node',
        'qremove_all_sync_node',
        #
        'displayln_malloc_stats',
        #
        'get_all_caches_info',
        'clear_all_caches',
        'clear_mem_cache',
        #
        'get_eigen_type',
        ]

__all__ += [
        'Gobble',
        #
        'QFile',
        'open_qfile',
        'open_qfile_str',
        'QarFile',
        'open_qar',
        'open_qar_info',
        'get_qar_multi_vol_max_size',
        'set_qar_multi_vol_max_size',
        'clean_up_qfile_map',
        'show_all_qfile',
        'does_regular_file_exist_qar',
        'does_file_exist_qar',
        'qar_build_index',
        'qar_create',
        'qar_extract',
        'qcopy_file',
        'list_qar',
        #
        'qcat',
        'qcat_bytes',
        'qtouch',
        'qappend',
        'qload_datatable',
        'compute_crc32',
        #
        'qar_build_index_info',
        'qar_create_info',
        'qar_extract_info',
        'qcopy_file_info',
        'qtouch_info',
        'qappend_info',
        'check_all_files_crc32_info',
        #
        'does_regular_file_exist_qar_sync_node',
        'does_file_exist_qar_sync_node',
        'qar_create_sync_node',
        'qar_extract_sync_node',
        'qcopy_file_sync_node',
        'qcat_sync_node',
        'qcat_bytes_sync_node',
        'qload_datatable_sync_node',
        ]

# ----------------------------------

import ctypes
import sys
import os
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | os.RTLD_GLOBAL)

lib_path = os.path.join(os.path.dirname(__file__),
                        'lib/libqlat-utils.so')

if not os.path.isfile(lib_path):
    lib_path = os.path.join(os.path.dirname(__file__),
                            'lib/libqlat-utils.dylib')

assert os.path.isfile(lib_path)

ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)

from .timer import *
from .cutils import *
from .types import *
from .coordinate import *
from .lat_data import *
from .rng_state import *
from .qar import *
from .mat import *

sys.setdlopenflags(flags)
