import ctypes
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

from qlat_utils.cp import *
from cqlat import *
from .cp import *

sys.setdlopenflags(flags)
