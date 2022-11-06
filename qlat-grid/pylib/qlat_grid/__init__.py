import ctypes
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

import cqlat_grid as cgi

sys.setdlopenflags(flags)

from qlat_grid.init import *
from qlat_grid.prop import *

from qlat_grid.get_include_dir import *
