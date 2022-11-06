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

cpp = meson.get_compiler('cpp')

py_mod = import('python')
py3 = py_mod.find_installation('python3')
message(py3.path())
message(py3.get_install_dir())

omp = dependency('openmp')
fftw = dependency('fftw3')
fftwf = dependency('fftw3f')
fftw_mpi = cpp.find_library('fftw3_mpi')
fftwf_mpi = cpp.find_library('fftw3f_mpi')
zlib = dependency('zlib')
math = cpp.find_library('m')

message('Collecting Grid information using \'grid-config\':')
grid_prefix = run_command('grid-config', '--prefix', check: true).stdout().strip()
message('grid: prefix =', grid_prefix)
grid_cxxflags = run_command('grid-config', '--cxxflags', check: true).stdout().strip().split()
grid_ldflags = run_command('grid-config', '--ldflags', check: true).stdout().strip().split()
grid_libs = run_command('grid-config', '--libs', check: true).stdout().strip().split()
message('Grid CXXFLAGS:', grid_cxxflags)
message('Grid LDFLAGS:', grid_ldflags)
message('Grid LIBS:', grid_libs)

grid = declare_dependency(
  compile_args: grid_cxxflags,
  link_args: grid_ldflags + grid_libs,
  dependencies: cpp.find_library('Grid'),
  ).as_system()

deps = [ omp, py3.dependency(), grid, fftw_mpi, fftwf_mpi, fftw, fftwf, zlib, math, ]

subdir('include')

subdir('pylib')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls bin/*', check: true)
scripts = c.stdout().strip().split('\n')
install_data(scripts, install_dir: 'bin', install_mode: 'rwxr-xr-x')
