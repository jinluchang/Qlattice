project('qlat-utils', 'cpp',
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
zlib = dependency('zlib')

deps = [ omp, py3.dependency(), zlib ]

subdir('include/qlat-utils')

subdir('pylib/cqlat_utils')

incdir = include_directories('include')

py3.extension_module('cqlat_utils',
  cpp_sources + config,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  # gnu_symbol_visibility: 'default',
  )

subdir('pylib/qlat_utils')

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls bin/*', check: true)
scripts = files(c.stdout().strip().split('\n'))
install_data(scripts, install_dir: get_option('bindir'), install_mode: 'rwxr-xr-x')

pkg = import('pkgconfig')
pkg.generate(
  name: 'qlat-utils',
  description: 'Qlattice utilities.',
  libraries: [ '-lz', '-lm' ],
  install_dir: get_option('libdir') / 'pkgconfig',
  )
