import sys

usage_message = """
Usage:
    python -m qlat_utils_config [--cxxflags] [--ldflags] [--libs] [--LD_LIBRARY_PATH] [--eigen-type]
""".strip()

if len(sys.argv) < 2:
    print(usage_message)
    exit()

from qlat_utils_config import get_include_list, get_lib_list, get_new_ld_library_path, get_eigen_type

output_args = []

for arg in sys.argv[1:]:
    if arg == "--cxxflags":
        output_args += [ f"-I{path}" for path in get_include_list() ]
    elif arg == "--ldflags":
        output_args += [ f"-L{path}" for path in get_lib_list() ]
    elif arg == "--libs":
        output_args += [ "-lqlat-utils", ]
    elif arg == "--LD_LIBRARY_PATH":
        output_args = [ get_new_ld_library_path(), ]
    elif arg == "--eigen-type":
        output_args = [ get_eigen_type(), ]

print(" ".join(output_args))

exit()
