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

if cpp.check_header('Eigen/Eigen', dependencies: deps)
  conf_data.set('QLAT_USE_GRID_EIGEN', false)
elif cpp.check_header('Grid/Eigen/Eigen', dependencies: deps)
  conf_data.set('QLAT_USE_GRID_EIGEN', true)
endif

config_file = configure_file(output: 'config.h', configuration: conf_data)
py3.install_sources(config_file, subdir: 'qlat_utils/include/qlat-utils')

header_files = files(
  'array.h',
  'assert.h',
  'cache.h',
  'compatible-endian.h',
  'compatible-omp.h',
  'complex.h',
  'coordinate-d.h',
  'coordinate.h',
  'core.h',
  'crc32.h',
  'eigen.h',
  'env.h',
  'handle.h',
  'lat-io.h',
  'mat.h',
  'matrix.h',
  'matrix-hmc.h',
  'mat-vec.h',
  'mvector.h',
  'py_convert.h',
  'qacc-func.h',
  'qacc.h',
  'qar-cache.h',
  'qar.h',
  'qcd-setting.h',
  'qutils-extra.h',
  'qutils.h',
  'qutils-io.h',
  'qutils-vec.h',
  'rng-state.h',
  'sha256.h',
  'show.h',
  'sprng-sha256.h',
  'timer.h',
  'types.h',
  'vector.h',
  'version.h',
  )
py3.install_sources(header_files, subdir: 'qlat_utils/include/qlat-utils')

all_header_files = header_files + config_file
