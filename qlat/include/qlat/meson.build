conf_data = configuration_data()

if not cpp.has_function('mallopt', prefix: '#include <malloc.h>')
  conf_data.set('QLAT_NO_MALLOPT', true)
endif

if not cpp.has_function('sysinfo', prefix: '#include <sys/sysinfo.h>')
  conf_data.set('QLAT_NO_SYSINFO', true)
endif

config_file = configure_file(output: 'config.h', configuration: conf_data)

py3.install_sources(config_file, subdir: 'qlat/include/qlat')

header_fns = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h | grep -v "^config.h$"',
  check: true).stdout().strip().split('\n')
header_files = files(header_fns)
py3.install_sources(header_files, subdir: 'qlat/include/qlat')

py3.install_sources('qlat-setup.h', subdir: 'qlat/include')

subdir('vector_utils')
subdir('lib')

all_header_files = header_files + config_file + vector_utils_header_files + lib_header_files
