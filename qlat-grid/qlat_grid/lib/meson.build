cppsources = files(
  'init.cpp',
  'prop.cpp',
  )
py3.install_sources(cppsources, subdir: 'qlat_grid/lib')

clib = shared_library('qlat-grid',
  cppsources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat_grid/lib',
  install_rpath: ':'.join(qlat_lib),
  )
