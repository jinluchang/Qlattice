c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py *.pxd', check: true)
pysources = files(c.stdout().strip().split('\n'))
py3.install_sources(pysources, subdir: 'qlat')

copy = []
foreach f : pysources
  copy += fs.copyfile(f)
endforeach

cplib = py3.extension_module('cp',
  files('cp.pyx'),
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat',
  )
