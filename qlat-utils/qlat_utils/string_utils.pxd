# Explicit UTF-8 conversion helpers.  These replace the former module-level
# implicit C-string conversions with visible calls.

from libcpp.string cimport string as std_string

cdef inline bytes as_bytes(object x):
    if isinstance(x, str):
        return (<str>x).encode('utf-8')
    if isinstance(x, bytes):
        return x
    raise TypeError("expected str or bytes, got " + type(x).__name__)

cdef inline str as_str_std(const std_string& s):
    return s.decode('utf-8')

cdef inline str as_str_c(const char* s):
    return (<bytes>s).decode('utf-8')
