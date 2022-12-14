conf_data = configuration_data()

if cpp.has_function('malloc_stats', prefix: '#include <malloc.h>', dependencies: deps)
  conf_data.set('QLAT_USE_MALLOC_STATS', true)
endif

if not cpp.has_function('aligned_alloc', prefix: '#include <cstdlib>', dependencies: deps)
  conf_data.set('QLAT_NO_ALIGNED_ALLOC', true)
endif

if cpp.check_header('endian.h', dependencies: deps)
  # do nothing
elif cpp.check_header('machine/endian.h', dependencies: deps)
  conf_data.set('QLAT_USE_MACHINE_ENDIAN_H', true)
else
  error('Cannot find endian.h header file.')
endif

config = configure_file(output: 'config.h', configuration: conf_data)
py3.install_sources(config, subdir: 'qlat_utils/include/qlat-utils')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h | grep -v "^config.h$"', check: true)
cpp_headers = files(c.stdout().strip().split('\n'))
py3.install_sources(cpp_headers, subdir: 'qlat_utils/include/qlat-utils')

subdir('lib')
