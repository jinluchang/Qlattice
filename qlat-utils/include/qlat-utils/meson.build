conf_data = configuration_data()

if cpp.has_function('malloc_stats', prefix: '#include <malloc.h>')
  conf_data.set('QLAT_USE_MALLOC_STATS', true)
endif

if cpp.check_header('endian.h')
  # do nothing
elif cpp.check_header('machine/endian.h')
  conf_data.set('QLAT_USE_MACHINE_ENDIAN_H', true)
else
  error('Cannot find endian.h header file.')
endif

config = configure_file(output: 'config.h', configuration: conf_data)

install_headers(config, subdir: 'qlat-utils')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h | grep -v "^config.h$"', check: true)
cpp_headers = files(c.stdout().strip().split('\n'))
install_headers(cpp_headers, subdir: 'qlat-utils')
