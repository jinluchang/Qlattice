#!/usr/bin/env python3

import sys

if not sys.version_info.major == 3 and sys.version_info.minor >= 6:
    print(sys.argv)
    print("You are using not supported Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

import os

def set_env(env, val):
    return f'export {env}="$setenv_prefix"/"{val}":"${env}"'

prefix = os.getenv("prefix")

assert prefix[:-1] != '/'

bin_dir = os.path.join(prefix, "bin")
include_dir = os.path.join(prefix, "include")
lib_dir = os.path.join(prefix, "lib")
lib64_dir = os.path.join(prefix, "lib64")
lib_pkg_config_dir = os.path.join(prefix, "lib/pkgconfig")
lib64_pkg_config_dir = os.path.join(prefix, "lib64/pkgconfig")
lib_python_dir = os.path.join(prefix, "lib/python3")

setenv_fn = os.path.join(prefix, "setenv.sh")

l_init = []
if "--keep" in sys.argv and os.path.isfile(setenv_fn):
    with open(setenv_fn, mode = "r") as f:
        l_init = f.readlines()
        l_init = [ v.rstrip() for v in l_init ]

l = []

l.append("#!/bin/bash")

l.append("")

l.append(f'setenv_prefix="{prefix}"')

if l_init:
    l.append("")
    l.append("# -------------------------------------------------------------------")
    l += l_init
    l.append("# -------------------------------------------------------------------")

l.append("")

if os.path.isdir(bin_dir):
    l.append(set_env("PATH", "bin"))
    l.append("")

if os.path.isdir(include_dir):
    l.append(set_env("C_INCLUDE_PATH", "include"))
    l.append(set_env("CPLUS_INCLUDE_PATH", "include"))
    l.append("")

if os.path.isdir(lib_dir):
    l.append(set_env("LD_RUN_PATH", "lib"))
    l.append(set_env("LD_LIBRARY_PATH", "lib"))
    l.append(set_env("LIBRARY_PATH", "lib"))
    l.append("")

if os.path.isdir(lib64_dir):
    l.append(set_env("LD_RUN_PATH", "lib64"))
    l.append(set_env("LD_LIBRARY_PATH", "lib64"))
    l.append(set_env("LIBRARY_PATH", "lib64"))
    l.append("")

if os.path.isdir(lib_pkg_config_dir):
    l.append(set_env("PKG_CONFIG_PATH", "lib/pkgconfig"))
    l.append("")

if os.path.isdir(lib64_pkg_config_dir):
    l.append(set_env("PKG_CONFIG_PATH", "lib64/pkgconfig"))
    l.append("")

if os.path.isdir(lib_python_dir):
    l.append(set_env("PYTHONPATH", "lib/python3"))
    l.append("")

l.append(f'unset setenv_prefix')

organize_env_path = f"""
if which python3 >/dev/null 2>&1 && which organize-env-path.py >/dev/null 2>&1 ; then
    source <(organize-env-path.py)
fi
"""

l.append(organize_env_path)

with open(setenv_fn, mode = "w") as f:
    f.write("\n".join(l))
