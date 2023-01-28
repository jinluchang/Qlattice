#!/usr/bin/env python3

import sys

if not sys.version_info.major == 3 and sys.version_info.minor >= 6:
    print(sys.argv)
    print("You are using not supported Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

import os
import sysconfig

user_scripts_path = sysconfig.get_path('scripts', f'{os.name}_user')
print(user_scripts_path)
