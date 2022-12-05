c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py *.pxd', check: true)
pysources = files(c.stdout().strip().split('\n'))
py3.install_sources(pysources, subdir: 'qlat_grid')

copy = []
foreach f : pysources
  copy += fs.copyfile(f)
endforeach

cplib = py3.extension_module('cp',
  files('cp.pyx'),
  dependencies: deps,
  include_directories: incdir,
  install: true,
  subdir: 'qlat_grid',
  # gnu_symbol_visibility: 'default',
  )

