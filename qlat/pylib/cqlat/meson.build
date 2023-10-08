# Need to update `exports.h` by running `bash update.sh` manually.

cpp_fns = run_command('bash', '-c', 'cd "$MESON_SOURCE_ROOT/$MESON_SUBDIR" ; ls *.cpp',
  check: true).stdout().strip().split('\n')
cpp_files = files(cpp_fns)

lib = py3.extension_module('cqlat',
  cpp_files + config_file,
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  install_rpath: ':'.join(qlat_utils_lib + [ '$ORIGIN/lib', ]),
  # gnu_symbol_visibility: 'default',
  )
