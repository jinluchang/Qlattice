#!/usr/bin/env python3

import os

def get_colon_list(str_colon_separated_list):
    l = str_colon_separated_list.split(':')
    if not l:
        return []
    elif l[-1] == "":
        l.pop()
    return l

def show_env(env):
    val = os.getenv(env)
    if val is None:
        return None
    l = get_colon_list(val)
    if not l:
        return f"{env} = []"
    l_str = "'\n    '".join(l)
    return f"{env} = [\n    '{l_str}'\n    ]"

def show_env_list(env_list):
    l = []
    for env in env_list:
        v = show_env(env)
        if v is not None:
            l.append(v)
    return "\n".join(l)

def show_env_var(env):
    val = os.getenv(env)
    if val is None:
        return None
    return f"{env} = '{val}'"

def show_env_var_list(env_list):
    l = []
    for env in env_list:
        v = show_env_var(env)
        if v is not None:
            l.append(v)
    return "\n".join(l)

all_env_list = [
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

all_env_var_list = [
        "prefix",
        "num_proc",
        "USE_COMPILER",
        "CC",
        "CXX",
        "MPICC",
        "MPICXX",
        "CFLAGS",
        "CXXFLAGS",
        "LDFLAGS",
        "LIBS",
        "LD_PRELOAD",
        ]

if __name__ == "__main__":
    print(show_env_var_list(all_env_var_list))
    print(show_env_list(all_env_list))
