project('qlat', 'cpp', 'cython',
  version: '0.83',
  license: 'GPL-3.0-or-later',
  default_options: [
    'warning_level=3',
    'cpp_std=c++17',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    'cython_language=cpp',
    ])

add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])

subdir('depend-qlat-utils')

cpp = meson.get_compiler('cpp')

py_mod = import('python')
py3 = py_mod.find_installation('python3', pure: false)
message(py3.full_path())
message(py3.get_install_dir())

fs = import('fs')

cxx = run_command('bash', '-c', 'echo "$CXX"', check: true).stdout().strip()
mpicxx = run_command('bash', '-c', 'echo "$MPICXX"', check: true).stdout().strip()

if cxx != '' and mpicxx == cxx
  message(f'cxx=\'@cxx@\' (use CXX compiler without additional MPI options.)')
  mpic = dependency('', required: false)
else
  message(f'cxx=\'@cxx@\' mpicxx=\'@mpicxx@\' (use meson\'s automatic MPI detection.)')
  mpic = dependency('mpi', language: 'cpp').as_system()
endif

fftw = dependency('fftw3').as_system()
fftwf = dependency('fftw3f').as_system()
message('fftw libdir', fftw.get_variable('libdir'))
message('fftwf libdir', fftwf.get_variable('libdir'))
fftw_all = [ fftw, fftwf, ]

cuba = cpp.find_library('cuba', required: false)
gsl = dependency('gsl').as_system()

deps = [ qlat_utils, mpic, fftw_all, cuba, gsl, ]

py3.install_sources('qlat_gpt.py')

subdir('qlat_config')
subdir('qlat')
subdir('cqlat')
subdir('auto_contractor')
subdir('qlat_scripts')
subdir('bin')
