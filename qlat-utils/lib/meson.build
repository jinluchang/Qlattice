c = run_command('bash', '-c', 'cd "$MESON_SOURCE_ROOT/$MESON_SUBDIR" ; ls *.cpp | grep -v version.cpp', check: true)
cppsources = files(c.stdout().strip().split('\n'))
py3.install_sources(cppsources, subdir: 'qlat_utils/lib')

version = custom_target('version',
  output: 'version.cpp',
  input: files([ 'version_gen.sh', 'version.cpp', ]),
  command: [ 'bash', '@INPUT@', '@OUTPUT@', ],
  install: false,
  build_always_stale: true,
  )

clib = static_library('qlat-utils',
  cppsources + version,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat_utils/lib',
  )
