project('qlat-ext', 'cpp', 'cython',
  version: '0.1',
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

py_mod = import('python')
py3 = py_mod.find_installation('python3')
message(py3.path())
message(py3.get_install_dir())

cxx = run_command('bash', '-c', 'echo "$CXX"', check: true).stdout().strip()
mpicxx = run_command('bash', '-c', 'echo "$MPICXX"', check: true).stdout().strip()

if cxx != '' and mpicxx == cxx
  message(f'cxx=\'@cxx@\' (use CXX compiler without additional MPI options.)')
  mpic = dependency('', required: false)
else
  message(f'cxx=\'@cxx@\' mpicxx=\'@mpicxx@\' (use meson\'s automatic MPI detection.)')
  mpic = dependency('mpi', language: 'cpp').as_system()
endif

omp = dependency('openmp').as_system()
zlib = dependency('zlib').as_system()

fftw = dependency('fftw3').as_system()
fftwf = dependency('fftw3f').as_system()
message('fftw libdir', fftw.get_variable('libdir'))
message('fftwf libdir', fftwf.get_variable('libdir'))
fftw_all = [ fftw, fftwf, ]

math = cpp.find_library('m')

numpy_include = run_command(py3.path(), '-c', 'import numpy as np ; print(np.get_include())', check: true).stdout().strip()
message('numpy include', numpy_include)

numpy = declare_dependency(
  include_directories:  include_directories(numpy_include),
  dependencies: [ py3.dependency(), ],
  )

if cpp.check_header('Eigen/Eigen')
  eigen = dependency('', required: false)
else
  eigen = dependency('eigen3').as_system()
endif

qlat_include = run_command(py3.path(), '-c', 'import qlat as q ; print("\\n".join(q.get_include_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat include', qlat_include)

qlat_lib = run_command(py3.path(), '-c', 'import qlat as q ; print("\\n".join(q.get_lib_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat lib', qlat_lib)

qlat_pxd = run_command(py3.path(), '-c', 'import qlat as q ; print("\\n".join(q.get_pxd_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat pxd', qlat_pxd)
qlat_pxd = files(qlat_pxd)

qlat_header = run_command(py3.path(), '-c', 'import qlat as q ; print("\\n".join(q.get_header_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat header', qlat_header)
qlat_header = files(qlat_header)

qlat = declare_dependency(
  include_directories: include_directories(qlat_include),
  dependencies: [
    py3.dependency().as_system(),
    cpp.find_library('qlat', dirs: qlat_lib),
    cpp.find_library('qlat-utils', dirs: qlat_lib),
    numpy, eigen, omp, fftw_all, mpic, zlib, math, ],
  )

deps = [ qlat, ]

subdir('include')
subdir('qlat_ext')

run_py = custom_target('run.py',
  output: 'run.py',
  input: files('run.py'),
  command: [ 'cp', '@INPUT@', '@OUTPUT@', ],
  install: true,
  install_dir: 'bin',
  )
