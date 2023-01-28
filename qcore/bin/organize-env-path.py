#!/usr/bin/env python3

import os

def organize_colon_list(str_colon_separated_list):
    if str_colon_separated_list.endswith(':'):
        str_colon_separated_list = str_colon_separated_list[:-1]
    l = str_colon_separated_list.split(':')
    l_new = []
    for v in l:
        if v not in l_new:
            l_new.append(v)
    return ':'.join(l_new)

def set_env(env):
    val = os.getenv(env)
    val = organize_colon_list(val)
    return f"export {env}='{val}'"

def set_env_list(env_list):
    l = []
    for env in env_list:
        l.append(set_env(env))
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

if __name__ == "__main__":
    print(set_env_list(all_env_list))
