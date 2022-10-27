project('qlat-utils', 'cpp',
  version: '0.2',
  license: 'GPL-3.0-or-later',
  default_options: [
    'warning_level=3',
    'cpp_std=c++14',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    ])

add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])

run_command('bash', '-c', 'cd $MESON_SOURCE_ROOT/pylib/cqlat_utils ; bash update.sh', check: true)

cpp = meson.get_compiler('cpp')

py_mod = import('python')
py3 = py_mod.find_installation('python3')
message(py3.path())
message(py3.get_install_dir())

omp = dependency('openmp')
zlib = dependency('zlib')

deps = [ omp, py3.dependency(), zlib ]

c = run_command('bash', '-c', 'ls $MESON_SOURCE_ROOT/include/qlat-utils/*.h', check: true)
headers = c.stdout().strip().split('\n')
install_headers(headers, subdir: 'qlat-utils')

c = run_command('bash', '-c', 'ls $MESON_SOURCE_ROOT/pylib/cqlat_utils/*.cpp', check: true)
sources = c.stdout().strip().split('\n')

incdir = [ include_directories('include') ]

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

py3.extension_module('cqlat_utils',
  sources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  # gnu_symbol_visibility: 'default',
  )

c = run_command('bash', '-c', 'ls $MESON_SOURCE_ROOT/pylib/qlat_utils/*.py', check: true)
pysources = c.stdout().strip().split('\n')
py3.install_sources(pysources, subdir: 'qlat_utils')

c = run_command('bash', '-c', 'ls $MESON_SOURCE_ROOT/bin/*', check: true)
scripts = c.stdout().strip().split('\n')
install_data(scripts, install_dir: 'bin')

pkg = import('pkgconfig')
pkg.generate(
  name: 'qlat-utils',
  description: 'Qlattice utilities.',
  libraries: [ '-lz', '-lm' ],
  )
