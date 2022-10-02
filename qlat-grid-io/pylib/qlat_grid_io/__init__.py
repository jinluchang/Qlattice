import ctypes
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

import cqlat_grid_io as cgi

sys.setdlopenflags(flags)

from qlat_grid_io.init import *
from qlat_grid_io.prop import *
