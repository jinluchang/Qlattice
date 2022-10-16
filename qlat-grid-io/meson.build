project('qlat-grid-io', 'cpp',
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

run_command('bash', '-c', 'cd "$MESON_SOURCE_ROOT"/pylib/cqlat_grid_io ; bash update.sh', check: true)

py_mod = import('python')
py3 = py_mod.find_installation('python3')
message(py3.path())
message(py3.get_install_dir())

omp = dependency('openmp')

qlat = dependency('qlat', method: 'pkg-config')

message('Collecting Grid information using meson option \'grid_prefix\':')
grid_prefix = get_option('grid_prefix')
grid_cxxflags = run_command(grid_prefix / 'bin' / 'grid-config', '--cxxflags', check: true).stdout().strip().split()
grid_ldflags = run_command(grid_prefix / 'bin' / 'grid-config', '--ldflags', check: true).stdout().strip().split()
grid_libs = run_command(grid_prefix / 'bin' / 'grid-config', '--libs', check: true).stdout().strip().split()
message('Grid CXXFLAGS:', grid_cxxflags)
message('Grid LDFLAGS:', grid_ldflags)
message('Grid LIBS:', grid_libs)

grid = declare_dependency(
  compile_args: grid_cxxflags,
  link_args: grid_ldflags + grid_libs + [ '-lGrid', ],
  )

deps = [ omp, py3.dependency(), grid, qlat ]

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/include/qlat-grid-io/*.h', check: true)
headers = c.stdout().strip().split('\n')
install_headers(headers, subdir: 'qlat-grid-io')

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/pylib/cqlat_grid_io/*.cpp', check: true)
sources = c.stdout().strip().split('\n')

incdir = include_directories('include')

lib = py3.extension_module('cqlat_grid_io',
  sources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  # gnu_symbol_visibility: 'default',
  )

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/pylib/qlat_grid_io/*.py', check: true)
pysources = c.stdout().strip().split('\n')
py3.install_sources(pysources, subdir: 'qlat_grid_io')

c = run_command('bash', '-c', 'ls $MESON_SOURCE_ROOT/bin/*', check: true)
scripts = c.stdout().strip().split('\n')
install_data(scripts, install_dir: 'bin')

pkg = import('pkgconfig')
pkg.generate(
  name: 'qlat-grid-io',
  description: 'IO utilities using Qlattice and Grid.',
  requires: [ 'qlat', ],
  )
