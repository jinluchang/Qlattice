import sys

usage_message = """
Usage:
    python3 -m qlat config [--cxxflags] [--ldflags] [--libs] [--LD_LIBRARY_PATH]
""".strip()

if len(sys.argv) < 2:
    print(usage_message)
    exit()

output_args = []

action = sys.argv[1]

if action == "config":
    from .get_include_dir import get_include_list, get_lib_list, get_new_ld_library_path
    for arg in sys.argv[2:]:
        if arg == "--cxxflags":
            output_args += [ f"-I{path}" for path in get_include_list() ]
        elif arg == "--ldflags":
            output_args += [ f"-L{path}" for path in get_lib_list() ]
        elif arg == "--libs":
            output_args += [ "-lqlat-utils", "-lqlat", ]
        elif arg == "--LD_LIBRARY_PATH":
            output_args = [ get_new_ld_library_path() ]
    print(" ".join(output_args))
    exit()

from .__init__ import *

print(usage_message)
exit()
