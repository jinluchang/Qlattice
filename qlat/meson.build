project('qlat', 'cpp',
  version: '0.1',
  license: 'GPL-3.0-or-later',
  default_options: [
    'warning_level=3',
    'cpp_std=c++14',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    ])

add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])

run_command('bash', '-c', 'cd "$MESON_SOURCE_ROOT"/pylib/cqlat ; bash update.sh', check: true)

cpp = meson.get_compiler('cpp')

py_mod = import('python')
py3 = py_mod.find_installation('python3')

message(py3.path())
message(py3.get_install_dir())

omp = dependency('openmp')

if get_option('use_cxx')
  message('use_cxx=true (use CXX compiler without additional MPI options.)')
  mpic = dependency('', required: false)
else
  message('use_cxx=false (use meson\'s automatic MPI detection.)')
  mpic = dependency('mpi', language: 'cpp')
endif

fftw = dependency('fftw3')
fftwf = dependency('fftw3f')
fftw_mpi = cpp.find_library('fftw3_mpi')
fftwf_mpi = cpp.find_library('fftw3f_mpi')
zlib = dependency('zlib')
math = cpp.find_library('m', required : false)

eigen = dependency('eigen3', required: false)

qutils = dependency('qlat-utils')

deps = [ omp, py3.dependency(), eigen, qutils, fftw_mpi, fftwf_mpi, fftw, fftwf, mpic, zlib, math, ]

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/include/qlat/*.h', check: true)
headers = c.stdout().strip().split('\n')
install_headers(headers, subdir: 'qlat')

install_headers('include/qlat/qlat-setup.h')

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/pylib/cqlat/*.cpp', check: true)
sources = c.stdout().strip().split('\n')

incdir = include_directories('include')

lib = py3.extension_module('cqlat',
  sources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  )

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/pylib/qlat/*.py', check: true)
pysources = c.stdout().strip().split('\n')
py3.install_sources(pysources, subdir: 'qlat')

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/pylib/auto_contractor/*.py', check: true)
pysources = c.stdout().strip().split('\n')
py3.install_sources(pysources, subdir: 'auto_contractor')

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/pylib/rbc_ukqcd_params/*.py', check: true)
pysources = c.stdout().strip().split('\n')
py3.install_sources(pysources, subdir: 'rbc_ukqcd_params')

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/pylib/*.py', check: true)
pysources = c.stdout().strip().split('\n')
py3.install_sources(pysources)

c = run_command('bash', '-c', 'ls $MESON_SOURCE_ROOT/bin/*', check: true)
scripts = c.stdout().strip().split('\n')
install_data(scripts, install_dir: 'bin')

pkg = import('pkgconfig')
pkg.generate(
  name: 'qlat',
  description: 'Qlattice - A simple lattice QCD library.',
  requires: [ 'fftw3', 'fftw3f', 'qlat-utils', ],
  libraries: [ '-lfftw3_mpi', '-lfftw3f_mpi', '-lz', '-lm' ],
  )
