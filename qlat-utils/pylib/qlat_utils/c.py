import ctypes
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

from cqlat_utils import *
from qlat_utils.cp import *

sys.setdlopenflags(flags)
