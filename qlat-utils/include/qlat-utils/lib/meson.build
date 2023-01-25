lib_header_fns = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h', check: true).stdout().strip().split('\n')
lib_header_files = files(lib_header_fns)
py3.install_sources(lib_header_files, subdir: 'qlat_utils/include/qlat-utils/lib')
