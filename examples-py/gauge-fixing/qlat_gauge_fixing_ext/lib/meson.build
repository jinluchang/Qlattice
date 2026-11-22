cppsources = files(
  'gauge-fixing.cpp',
  )
py3.install_sources(cppsources, subdir: 'qlat_gauge_fixing_ext/lib')

clib = shared_library('qlat-gauge-fixing',
  cppsources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat_gauge_fixing_ext/lib',
  install_rpath: ':'.join(qlat_lib),
  )
