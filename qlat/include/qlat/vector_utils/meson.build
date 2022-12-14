c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h', check: true)
vector_utils_cpp_headers = files(c.stdout().strip().split('\n'))
py3.install_sources(vector_utils_cpp_headers, subdir: 'qlat/include/qlat/vector_utils')
