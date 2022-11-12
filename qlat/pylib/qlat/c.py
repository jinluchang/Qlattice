import ctypes
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

from cqlat_utils import *
from cqlat import *

sys.setdlopenflags(flags)
