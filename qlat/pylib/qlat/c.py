from qlat_utils.c import *

import ctypes
import sys
import os
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

ctypes.cdll.LoadLibrary(
        os.path.join(
            os.path.dirname(__file__),
            'lib/libqlat.so'))

from cqlat import *
from .cp import *

sys.setdlopenflags(flags)
