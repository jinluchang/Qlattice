from qlat.c import *

import ctypes
import os
import sys

flags = sys.getdlopenflags()
sys.setdlopenflags(flags | os.RTLD_GLOBAL)

lib_path = os.path.join(os.path.dirname(__file__), "lib/libqlat-gauge-fixing.so")

if not os.path.isfile(lib_path):
    lib_path = os.path.join(os.path.dirname(__file__), "lib/libqlat-gauge-fixing.dylib")

assert os.path.isfile(lib_path)

ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)

from .gauge_fixing import *

sys.setdlopenflags(flags)
