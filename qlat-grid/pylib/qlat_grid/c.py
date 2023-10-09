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

ctypes.cdll.LoadLibrary(lib_path)

from .cp import *

sys.setdlopenflags(flags)
