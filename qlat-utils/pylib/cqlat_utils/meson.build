# Need to update `exports.h` by running `bash update.sh` manually.

cpp_fns = run_command(
  'bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.cpp',
  check: true,
  ).stdout().strip().split('\n')
cpp_files = files(cpp_fns)

py3.extension_module('cqlat_utils',
  cpp_files + config_file,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  # gnu_symbol_visibility: 'default',
  )
