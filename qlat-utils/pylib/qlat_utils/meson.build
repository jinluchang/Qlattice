c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py *.pxd', check: true)
qlat_utils_pysources = files(c.stdout().strip().split('\n'))
py3.install_sources(qlat_utils_pysources, subdir: 'qlat_utils')

copy = []
foreach f : qlat_utils_pysources
  copy += fs.copyfile(f)
endforeach

cplib = py3.extension_module('cp',
  files('cp.pyx'),
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat_utils',
  )
