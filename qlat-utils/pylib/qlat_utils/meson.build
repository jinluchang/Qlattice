c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py', check: true)
qlat_utils_pysources = files(c.stdout().strip().split('\n'))
py3.install_sources(qlat_utils_pysources, subdir: 'qlat_utils')
