c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h', check: true)
libhsources = files(c.stdout().strip().split('\n'))
py3.install_sources(libhsources, subdir: 'qlat/include/qlat/lib')
