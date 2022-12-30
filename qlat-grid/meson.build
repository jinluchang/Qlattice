project('qlat-grid', 'cpp', 'cython',
  version: '0.18',
  license: 'GPL-3.0-or-later',
  default_options: [
    'warning_level=3',
    'cpp_std=c++14',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    'cython_language=cpp',
    ])

add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])

cpp = meson.get_compiler('cpp')

fs = import('fs')

py_mod = import('python')
py3 = py_mod.find_installation('python3')
message(py3.path())
message(py3.get_install_dir())

message('Collecting Grid information using \'grid-config\':')
grid_prefix = run_command('grid-config', '--prefix', check: true).stdout().strip()
message('Grid: prefix =', grid_prefix)
message('Grid: CXX =', run_command('grid-config', '--cxx', check: true).stdout().strip())
message('Grid: CXXLD =', run_command('grid-config', '--cxxld', check: true).stdout().strip())
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
    v = f.substring(2)
    if v not in grid_include_directories
      grid_include_directories += v
    endif
  else
    if f not in grid_compile_args
      grid_compile_args += f
    endif
  endif
endforeach

qlat_lib_names = [ 'z', 'm', 'fftw3', 'fftw3f', ]

grid_link_args = []
grid_lib_dirs = []
grid_lib_names = [ 'Grid', ]
foreach f : grid_ldflags + grid_libs
  if f.startswith('-L')
    v = f.substring(2)
    if v not in grid_lib_dirs
      grid_lib_dirs += v
    endif
  elif f.startswith('-l')
    v = f.substring(2)
    if v not in grid_lib_names + qlat_lib_names
      grid_lib_names += v
    endif
  else
    if f not in grid_link_args
      grid_link_args += f
    endif
  endif
endforeach

message('Grid processed info:')
message('Grid: compile_args =', grid_compile_args)
message('Grid: include_directories =', grid_include_directories)
message('Grid: link_args =', grid_link_args)
message('Grid: lib_dirs =', grid_lib_dirs)
message('Grid: lib_names =', grid_lib_names)

grid_deps = []
grid_link_args_not_found_libs = []
foreach lib : grid_lib_names
  dep = cpp.find_library(lib, dirs: grid_lib_dirs, required: false)
  if dep.found()
    grid_deps += dep
  else
    grid_link_args_not_found_libs += '-l' + lib
  endif
endforeach

message('Grid: link_args_not_found_libs =', grid_link_args_not_found_libs)

grid = declare_dependency(
  include_directories: grid_include_directories,
  compile_args: grid_compile_args,
  link_args: grid_link_args + grid_link_args_not_found_libs,
  dependencies: grid_deps,
  ).as_system()

omp = dependency('openmp').as_system()
zlib = dependency('zlib').as_system()

fftw = dependency('fftw3').as_system()
fftwf = dependency('fftw3f').as_system()
message('fftw libdir', fftw.get_variable('libdir'))
message('fftwf libdir', fftwf.get_variable('libdir'))
fftw_all = [ fftw, fftwf, ]

# fftw_mpi = cpp.find_library('fftw3_mpi', dirs: fftw.get_variable('libdir'))
# fftwf_mpi = cpp.find_library('fftw3f_mpi', dirs: fftwf.get_variable('libdir'))
# fftw_all += [ fftw_mpi, fftwf_mpi, ]

math = cpp.find_library('m')

numpy_include = run_command(py3.path(), '-c', 'import numpy as np ; print(np.get_include())', check: true).stdout().strip()
message('numpy include', numpy_include)

numpy = declare_dependency(
  include_directories:  include_directories(numpy_include),
  dependencies: [ py3.dependency(), ],
  ).as_system()

if cpp.check_header('Grid/Eigen/Eigen', dependencies: grid)
  eigen = dependency('', required: false)
else
  error('Cannot find Eigen headers provided by Grid.')
endif

qlat_include = run_command(py3.path(), '-c', 'import qlat as q ; print("\\n".join(q.get_include_list()))', env: environment({'q_verbose': '-1'}), check: true).stdout().strip().split('\n')
message('qlat include', qlat_include)

qlat_lib = run_command(py3.path(), '-c', 'import qlat as q ; print("\\n".join(q.get_lib_list()))', env: environment({'q_verbose': '-1'}), check: true).stdout().strip().split('\n')
message('qlat lib', qlat_lib)

qlat = declare_dependency(
  include_directories: include_directories(qlat_include),
  dependencies: [
    py3.dependency().as_system(),
    cpp.find_library('qlat', dirs: qlat_lib),
    cpp.find_library('qlat-utils', dirs: qlat_lib),
    numpy, eigen, omp, fftw_all, zlib, math, ],
  )

deps = [ qlat, grid, ]

subdir('include')

subdir('pylib')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls bin/*', check: true)
scripts = files(c.stdout().strip().split('\n'))
install_data(scripts, install_dir: get_option('bindir'), install_mode: 'rwxr-xr-x')
