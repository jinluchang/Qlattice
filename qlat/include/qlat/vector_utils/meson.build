vector_utils_header_fns = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h',
  check: true).stdout().strip().split('\n')
vector_utils_header_files = files(vector_utils_header_fns)
py3.install_sources(vector_utils_header_files, subdir: 'qlat/include/qlat/vector_utils')
