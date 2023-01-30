#!/usr/bin/env python3

# Usage:
# find-library.py lib-name-1 lib-name-2 ...
# return first prefix of the path in LIBRARY_PATH (endswith /lib or /lib64) that contains one of the lib-names

import sys

if not (sys.version_info.major == 3 and sys.version_info.minor >= 6):
    print(sys.argv)
    print("You are using not supported Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

import os
import glob

def get_colon_list(str_colon_separated_list):
    l = str_colon_separated_list.split(':')
    if not l:
        return []
    elif l[-1] == "":
        l.pop()
    return l

def find_library(*lib_list):
    """find a first path that contains one of lib_list"""
    val = os.getenv("LIBRARY_PATH")
    if val is None:
        return ""
    l = get_colon_list(val)
    for p in l:
        for lib in lib_list:
            if glob.glob(f"{p}/{lib}*"):
                if p.endswith("/lib"):
                    return p[:-len('/lib')]
                elif p.endswith("/lib64"):
                    return p[:-len('/lib64')]
    return ""

if __name__ == "__main__":
    print(find_library(*sys.argv[1:]))
