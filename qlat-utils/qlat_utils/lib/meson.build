version = custom_target('version',
  output: 'version.cpp',
  input: files([ 'version_gen.sh', 'version.cpp', ]),
  command: [ 'bash', '@INPUT@', '@OUTPUT@', ],
  install: false,
  build_always_stale: true,
  )
cppsources = files(
  'env.cpp',
  'timer.cpp',
  'lat-io.cpp',
  'mat.cpp',
  'qar.cpp',
  'qutils-io.cpp',
  'rng-state.cpp',
  )
py3.install_sources(cppsources, subdir: 'qlat_utils/lib')

clib = shared_library('qlat-utils',
  cppsources + version,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat_utils/lib',
  )
