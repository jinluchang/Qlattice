project('qlat-ext', 'cpp',
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

run_command('bash', '-c', 'cd "$MESON_SOURCE_ROOT"/pylib/cqlat_ext ; bash update.sh', check: true)

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

qlat = dependency('qlat', method: 'pkg-config')

deps = [ omp, py3.dependency(), mpic, qlat ]

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/pylib/cqlat_ext/*.cpp', check: true)
sources = c.stdout().strip().split('\n')

incdir = include_directories('include')

lib = py3.extension_module('cqlat_ext',
  sources,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  # gnu_symbol_visibility: 'default',
  )

pkg = import('pkgconfig')
pkg.generate(
  name: 'qlat-ext',
  description: 'Qlattice - A simple lattice QCD library.',
  requires: [ 'qlat', ],
  )
