cppsources = files(
  'init.cpp',
  'prop.cpp',
  )
py3.install_sources(cppsources, subdir: 'qlat_cps/lib')

clib = shared_library('qlat-cps',
  cppsources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat_cps/lib',
  install_rpath: ':'.join(qlat_lib),
  )
