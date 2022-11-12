project('qlat-grid', 'cpp',
  version: '0.5',
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

message('Collecting Grid information using \'grid-config\':')
grid_prefix = run_command('grid-config', '--prefix', check: true).stdout().strip()
message('Grid: prefix =', grid_prefix)
grid_cxxflags = run_command('grid-config', '--cxxflags', check: true).stdout().strip().split()
grid_ldflags = run_command('grid-config', '--ldflags', check: true).stdout().strip().split()
grid_libs = run_command('grid-config', '--libs', check: true).stdout().strip().split()
message('Grid: CXXFLAGS =', grid_cxxflags)
message('Grid: LDFLAGS =', grid_ldflags)
message('Grid: LIBS =', grid_libs)

grid_compile_args = []
grid_include_directories = []
foreach f : grid_cxxflags
  if f.startswith('-I')
    grid_include_directories += f.substring(2)
  else
    grid_compile_args += f
  endif
endforeach

grid_link_args = []
grid_lib_dirs = []
grid_lib_names = [ 'Grid', ]
foreach f : grid_ldflags + grid_libs
  if f.startswith('-L')
    grid_lib_dirs += f.substring(2)
  elif f.startswith('-l')
    grid_lib_names += f.substring(2)
  else
    grid_link_args += f
  endif
endforeach

message('Grid processed info:')
message('Grid: compile_args =', grid_compile_args)
message('Grid: include_directories =', grid_include_directories)
message('Grid: link_args =', grid_link_args)
message('Grid: lib_dirs =', grid_lib_dirs)
message('Grid: lib_names =', grid_lib_names)

grid_deps = []
foreach lib : grid_lib_names
  grid_deps += cpp.find_library(lib, dirs: grid_lib_dirs)
endforeach

grid = declare_dependency(
  include_directories: grid_include_directories,
  compile_args: grid_compile_args,
  link_args: grid_link_args,
  dependencies: grid_deps,
  ).as_system()

omp = dependency('openmp').as_system()
zlib = dependency('zlib').as_system()

fftw = dependency('fftw3').as_system()
fftwf = dependency('fftw3f').as_system()
message('fftw libdir', fftw.get_variable('libdir'))
message('fftwf libdir', fftwf.get_variable('libdir'))
fftw_mpi = cpp.find_library('fftw3_mpi', dirs: fftw.get_variable('libdir'))
fftwf_mpi = cpp.find_library('fftw3f_mpi', dirs: fftwf.get_variable('libdir'))

math = cpp.find_library('m')

qlat_include = run_command('qlat-include', check: true).stdout().strip().split('\n')
message('qlat include', qlat_include)

qlat = declare_dependency(
  include_directories: include_directories(qlat_include),
  dependencies: [ py3.dependency().as_system(), omp, fftw_mpi, fftwf_mpi, fftw, fftwf, zlib, math, ],
  )

deps = [ qlat, grid, ]

if not cpp.check_header('Eigen/Eigen')
  eigen = dependency('eigen3').as_system()
  deps += [ eigen, ]
endif

subdir('include')

subdir('pylib')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls bin/*', check: true)
scripts = files(c.stdout().strip().split('\n'))
install_data(scripts, install_dir: get_option('bindir'), install_mode: 'rwxr-xr-x')
