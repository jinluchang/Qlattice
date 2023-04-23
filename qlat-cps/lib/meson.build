c = run_command('bash', '-c', 'cd "$MESON_SOURCE_ROOT/$MESON_SUBDIR" ; ls *.cpp', check: true)
cppsources = files(c.stdout().strip().split('\n'))
py3.install_sources(cppsources, subdir: 'qlat-cps/lib')

clib = static_library('qlat-cps',
  cppsources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat-cps/lib',
  )
