c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h | grep -v "^config.h$"', check: true)
headers = files(c.stdout().strip().split('\n'))
install_headers(headers, subdir: 'qlat-grid-io')
