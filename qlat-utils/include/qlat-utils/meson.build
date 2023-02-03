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

if cpp.check_header('Grid/Eigen/Eigen')
  conf_data.set('QLAT_USE_GRID_EIGEN', true)
endif

config_file = configure_file(output: 'config.h', configuration: conf_data)
py3.install_sources(config_file, subdir: 'qlat_utils/include/qlat-utils')

header_fns = run_command(
  'bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.h | grep -v "^config.h$"',
  check: true,
  ).stdout().strip().split('\n')
header_files = files(header_fns)
py3.install_sources(header_files, subdir: 'qlat_utils/include/qlat-utils')

subdir('lib')

all_header_files = header_files + config_file + lib_header_files
