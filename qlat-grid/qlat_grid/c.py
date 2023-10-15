from qlat.c import *

import ctypes
import sys
import os
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | os.RTLD_GLOBAL)

lib_path = os.path.join(os.path.dirname(__file__),
                        'lib/libqlat-grid.so')

if not os.path.isfile(lib_path):
    lib_path = os.path.join(os.path.dirname(__file__),
                            'lib/libqlat-grid.dylib')

assert os.path.isfile(lib_path)

ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)

from .cp import *

sys.setdlopenflags(flags)
