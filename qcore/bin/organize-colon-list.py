#!/usr/bin/env python3

import sys

if not (sys.version_info.major == 3 and sys.version_info.minor >= 6):
    print(sys.argv)
    print("You are using not supported Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

def organize_colon_list(str_colon_separated_list):
    l = str_colon_separated_list.split(':')
    if not l:
        return []
    elif l[-1] == "":
        l.pop()
    l_new = []
    for v in l:
        if v not in l_new:
            l_new.append(v)
    return ':'.join(l_new)

program_name = sys.argv[0]

usage = f"""
Usage:\n
{program_name} colon-separated-list\n
Print a reorganized colon separated list to output.
"""

if __name__ == "__main__":
    #
    if len(sys.argv) != 2:
        print(usage)
        exit()
    #
    str_colon_separated_list = sys.argv[1]
    print(organize_colon_list(str_colon_separated_list))
