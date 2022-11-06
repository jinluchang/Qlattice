c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py', check: true)
pysources = files(c.stdout().strip().split('\n'))
py3.install_sources(pysources, subdir: 'rbc_ukqcd_params')
