cppsources = files(
  'core.cpp',
  'env.cpp',
  'field.cpp',
  'field-expand.cpp',
  'field-fft.cpp',
  'field-shuffle.cpp',
  'fields-io.cpp',
  'mpi.cpp',
  'selected-field.cpp',
  'selected-points.cpp',
  'utils-io.cpp',
  )
py3.install_sources(cppsources, subdir: 'qlat/lib')

clib = shared_library('qlat',
  cppsources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat/lib',
  install_rpath: ':'.join(qlat_utils_lib),
  )
