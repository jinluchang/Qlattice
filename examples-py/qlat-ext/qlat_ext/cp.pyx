from .everything cimport *
from qlat.everything cimport *

cdef void hello_world_c():
    cdef Timer timer = Timer(b"hello_world")
    timer.start()
    print("hello world")
    timer.stop()

def hello_world():
    hello_world_c()

def hello_cpy(x):
    print(f"hello {x}")
