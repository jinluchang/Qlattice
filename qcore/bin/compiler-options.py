#!/usr/bin/env python3

import sys

if not (sys.version_info.major == 3 and sys.version_info.minor >= 6):
    print(sys.argv)
    print("You are using not supported Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

import subprocess as p

def process_remove_arg(argv):
    name = "--wrapper-remove-arg="
    name_len = len(name)
    args_to_remove = []
    for arg in argv:
        if arg[:name_len] == name:
            args_to_remove.append(arg)
            args_to_remove.append(arg[name_len:])
    for x in args_to_remove:
        while x in argv:
            argv.remove(x)

argv = sys.argv.copy()
process_remove_arg(argv)
# status = p.run(argv[1:])
# sys.exit(status.returncode)
status = p.call(argv[1:])
sys.exit(status)
