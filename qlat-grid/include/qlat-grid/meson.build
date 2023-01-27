header_fns = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h | grep -v "^config.h$"',
  check: true).stdout().strip().split('\n')
header_files = files(header_fns)
py3.install_sources(header_files, subdir: 'qlat_grid/include/qlat-grid')

all_header_files = header_files
