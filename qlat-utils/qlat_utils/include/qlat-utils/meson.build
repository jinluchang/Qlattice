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

if eigen_type == 'grid'
  conf_data.set('QLAT_USE_GRID_EIGEN', true)
elif eigen_type == 'plain' or eigen_type == 'system'
  conf_data.set('QLAT_USE_GRID_EIGEN', false)
else
  error('eigen_type not correct')
endif

if quadmath.found()
  conf_data.set('__QLAT_NO_FLOAT128__', false)
else
  conf_data.set('__QLAT_NO_FLOAT128__', true)
endif

conf_data.set('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')

# conf_data.set('_GLIBCXX_DEBUG', true)
# conf_data.set('_GLIBCXX_DEBUG_PEDANTIC', true)
# conf_data.set('_GLIBCXX_DEBUG_BACKTRACE', true)

config_file = configure_file(output: 'config.h', configuration: conf_data)
py3.install_sources(config_file, subdir: 'qlat_utils/include/qlat-utils')

header_files = files(
  'array.h',
  'qassert.h',
  'cache.h',
  'compatible-endian.h',
  'compatible-omp.h',
  'complex.h',
  'coordinate-d.h',
  'coordinate.h',
  'core.h',
  'crc32.h',
  'eigen.h',
  'eigen-sys.h',
  'env.h',
  'qendian.h',
  'handle.h',
  'lat-io.h',
  'mat.h',
  'matrix.h',
  'matrix-hmc.h',
  'mat-vec.h',
  'mpi-auto.h',
  'mvector.h',
  'py_convert.h',
  'qacc-func.h',
  'qacc.h',
  'qar.h',
  'qcd-setting.h',
  'utils-extra.h',
  'utils.h',
  'utils-io.h',
  'utils-vec.h',
  'rng-state.h',
  'sha256.h',
  'show.h',
  'sprng-sha256.h',
  'timer.h',
  'types.h',
  'vector.h',
  'version.h',
  'warning.h',
  'qacc-translator.h',
  'utils-ldouble.h',
  )
py3.install_sources(header_files, subdir: 'qlat_utils/include/qlat-utils')

all_header_files = header_files + config_file
