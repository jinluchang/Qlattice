# Need to update `exports.h` by running `bash update.sh` manually.

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.cpp', check: true)
cpp_sources = files(c.stdout().strip().split('\n'))

pylib = py3.extension_module('cqlat_ext',
  cpp_sources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  # gnu_symbol_visibility: 'default',
  )
