# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT
cimport numpy

import numpy as np
import functools

include "elem_type.inc.pyx"

include "buffer.inc.pyx"

include "timer.inc.pyx"

include "mat.inc.pyx"
