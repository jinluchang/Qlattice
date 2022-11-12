project('qlat', 'cpp',
  version: '0.6',
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

if get_option('use_cxx')
  message('use_cxx=true (use CXX compiler without additional MPI options.)')
  mpic = dependency('', required: false)
else
  message('use_cxx=false (use meson\'s automatic MPI detection.)')
  mpic = dependency('mpi', language: 'cpp').as_system()
endif

omp = dependency('openmp').as_system()
zlib = dependency('zlib').as_system()

fftw = dependency('fftw3').as_system()
fftwf = dependency('fftw3f').as_system()
message('fftw libdir', fftw.get_variable('libdir'))
message('fftwf libdir', fftwf.get_variable('libdir'))
fftw_mpi = cpp.find_library('fftw3_mpi', dirs: fftw.get_variable('libdir'))
fftwf_mpi = cpp.find_library('fftw3f_mpi', dirs: fftwf.get_variable('libdir'))

math = cpp.find_library('m')

qlat_utils_include = run_command(py3.path(), '-c', 'import qlat_utils as q ; print("\\n".join(q.get_include_list()))', check: true).stdout().strip().split('\n')
message('qlat_utils include', qlat_utils_include)

qlat_utils = declare_dependency(
  include_directories:  include_directories(qlat_utils_include),
  dependencies: [ py3.dependency(), omp, zlib, math, ],
  )

deps = [ qlat_utils, fftw_mpi, fftwf_mpi, fftw, fftwf, mpic, ]

if not cpp.check_header('Eigen/Eigen')
  eigen = dependency('eigen3').as_system()
  deps += [ eigen, ]
endif

subdir('include')

subdir('pylib')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls bin/*', check: true)
scripts = files(c.stdout().strip().split('\n'))
install_data(scripts, install_dir: get_option('bindir'), install_mode: 'rwxr-xr-x')
