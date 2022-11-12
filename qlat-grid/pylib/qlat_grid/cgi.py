import ctypes
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

import cqlat_grid as cgi

sys.setdlopenflags(flags)
