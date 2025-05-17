# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat.all cimport *
from . cimport everything as cc

cdef void hello_world_c():
    cdef cc.Timer timer = cc.Timer(b"hello_world")
    timer.start()
    print("hello world")
    timer.stop()

def hello_world():
    hello_world_c()

def hello_cpy(x):
    print(f"hello {x}")

def hello():
    return cc.hello()
