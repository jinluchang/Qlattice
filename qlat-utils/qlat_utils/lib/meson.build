cppsource_fns = [
  'lat-io.cpp',
  'mat.cpp',
  'qar.cpp',
  'rng-state.cpp',
  'utils-io.cpp',
  ]
cppsources = files(cppsource_fns)
py3.install_sources(cppsources, subdir: 'qlat_utils/lib')

version = custom_target('version',
  output: 'version.cpp',
  input: files([ 'version_gen.sh', 'version.cpp', ]),
  command: [ 'bash', '@INPUT@', '@OUTPUT@', ],
  install: false,
  build_always_stale: true,
  )

clib = shared_library('qlat-utils',
  cppsources + version,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat_utils/lib',
  )
