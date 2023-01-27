#!/usr/bin/env python3

import sys

def organize(str_colon_separated_list):
    if str_colon_separated_list.endswith(':'):
        str_colon_separated_list = str_colon_separated_list[:-1]
    l = str_colon_separated_list.split(':')
    l_new = []
    for v in l:
        if v not in l_new:
            l_new.append(v)
    return ':'.join(l_new)

argv = sys.argv
program_name = argv[0]

usage = f"""
Usage:

{program_name} colon-separated-list

Print a reorganized colon separated list to output.
"""

if len(argv) != 2:
    print(usage)
    exit()

str_colon_separated_list = argv[1]
print(organize(str_colon_separated_list))
