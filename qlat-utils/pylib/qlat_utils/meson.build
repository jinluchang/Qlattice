c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py *.pxd', check: true)
qlat_utils_pysources = files(c.stdout().strip().split('\n'))
py3.install_sources(qlat_utils_pysources, subdir: 'qlat_utils')

copy = []
foreach f : qlat_utils_pysources
  copy += fs.copyfile(f)
endforeach

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.cpp', check: true)
cppsources = files(c.stdout().strip().split('\n'))

cplib = py3.extension_module('cp',
  files('cp.pyx') + cppsources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  subdir: 'qlat_utils',
  # gnu_symbol_visibility: 'default',
  )
