#!/usr/bin/env python3

import sys

if not (sys.version_info.major == 3 and sys.version_info.minor >= 6):
    print(sys.argv)
    print("You are using not supported Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

import os
import sysconfig
import site

def get_colon_list(str_colon_separated_list):
    l = str_colon_separated_list.split(':')
    if not l:
        return []
    elif l[-1] == "":
        l.pop()
    return l

def unique_list(l):
    l_new = []
    for v in l:
        if v not in l_new:
            l_new.append(v)
    return l_new

def organize_colon_list(str_colon_separated_list):
    return ':'.join(unique_list(get_colon_list(str_colon_separated_list)))

def set_env(env):
    val = os.getenv(env)
    if val is None:
        return None
    if env == "PATH":
        l = get_colon_list(val)
        user_scripts_path = sysconfig.get_path('scripts', f'{os.name}_user')
        if user_scripts_path not in l:
            l = [ user_scripts_path, ] + l
        val = ':'.join(unique_list(l))
    elif env == "PYTHONPATH":
        l = get_colon_list(val)
        l_site = site.getsitepackages()
        l_new = []
        for v in l:
            if v not in l_site:
                l_new.append(v)
        val = ':'.join(unique_list(l_new))
    else:
        val = organize_colon_list(val)
    return f"export {env}='{val}'"

def set_env_list(env_list):
    l = []
    for env in env_list:
        v = set_env(env)
        if v is not None:
            l.append(v)
    return "\n".join(l)

all_env_list = [
        "SETENV_PATH",
        "PATH",
        "PYTHONPATH",
        "LD_RUN_PATH",
        "LD_LIBRARY_PATH",
        "LIBRARY_PATH",
        "CPATH",
        "C_INCLUDE_PATH",
        "CPLUS_INCLUDE_PATH",
        "PKG_CONFIG_PATH",
        ]

if __name__ == "__main__":
    print(set_env_list(all_env_list))
