c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h | grep -v "^config.h$"', check: true)
headers = files(c.stdout().strip().split('\n'))
py3.install_sources(headers, subdir: 'qlat_grid/include/qlat-grid')
