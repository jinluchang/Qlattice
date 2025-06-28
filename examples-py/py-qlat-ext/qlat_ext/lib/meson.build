cppsources = files(
  'cp.cpp',
  )
py3.install_sources(cppsources, subdir: 'qlat_ext/lib')

clib = shared_library('qlat-ext',
  cppsources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat_ext/lib',
  install_rpath: ':'.join(qlat_lib),
  )
