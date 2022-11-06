project('qlat-cpp', 'cpp',
  version: '0.1',
  license: 'GPL-3.0-or-later',
  default_options: [
    'warning_level=0',
    'cpp_std=c++14',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    ])

add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])

cpp = meson.get_compiler('cpp')

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

qlat_include = run_command('qlat-include', check: true).stdout().strip().split('\n')
message('qlat include', qlat_include)

qlat = declare_dependency(
  include_directories: include_directories(qlat_include),
  dependencies: [ omp, fftw_mpi, fftwf_mpi, fftw, fftwf, mpic, zlib, math, ],
  )

deps = [ qlat, ]

if not cpp.check_header('Eigen/Eigen')
  eigen = dependency('eigen3').as_system()
  deps += [ eigen, ]
endif

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.cpp', check: true)
cpp_sources = c.stdout().strip().split('\n')

executable('qlat.x',
  cpp_sources,
  dependencies: deps,
  install: true,
  )

run_target('run',
  command: 'run.sh')
