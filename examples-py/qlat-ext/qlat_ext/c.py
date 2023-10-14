import ctypes
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

from .hello import *

sys.setdlopenflags(flags)
