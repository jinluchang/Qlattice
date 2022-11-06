conf_data = configuration_data()

if not cpp.has_function('mallopt', prefix: '#include <malloc.h>')
  conf_data.set('QLAT_NO_MALLOPT', true)
endif

if not cpp.has_function('sysinfo', prefix: '#include <sys/sysinfo.h>')
  conf_data.set('QLAT_NO_SYSINFO', true)
endif

config = configure_file(output: 'config.h', configuration: conf_data)

install_headers(config, subdir: 'qlat')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h | grep -v "^config.h$"', check: true)
cpp_headers = files(c.stdout().strip().split('\n'))
install_headers(cpp_headers, subdir: 'qlat')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls vector_utils/*.h', check: true)
vector_utils_cpp_headers = files(c.stdout().strip().split('\n'))
install_headers(vector_utils_cpp_headers, subdir: 'qlat/vector_utils')

install_headers('qlat-setup.h')
