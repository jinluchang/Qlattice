from qlat_utils.c import *

import ctypes
import sys
import os
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

lib_path = os.path.join(os.path.dirname(__file__),
                        'lib/libqlat.so')

if not os.path.isfile(lib_path):
    lib_path = os.path.join(os.path.dirname(__file__),
                            'lib/libqlat.dylib')

assert os.path.isfile(lib_path)

ctypes.cdll.LoadLibrary(lib_path)

from cqlat import *
from .cp import *

sys.setdlopenflags(flags)
