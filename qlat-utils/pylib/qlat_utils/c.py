import ctypes
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

from .cp import *
from .cpa import *

sys.setdlopenflags(flags)
