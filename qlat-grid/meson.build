project('qlat-grid', 'cpp', 'cython',
  version: '0.41',
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

subdir('depend-qlat')
subdir('depend-grid')

py_mod = import('python')
py3 = py_mod.find_installation('python3', pure: false)
message(py3.path())
message(py3.get_install_dir())

deps = [ qlat, grid, ]

subdir('qlat_grid')
