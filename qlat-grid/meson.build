project('qlat-grid', 'cpp', 'cython',
  version: '0.83',
  license: 'GPL-3.0-or-later',
  default_options: [
    # 'warning_level=0',
    'cpp_std=c++17',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    'cython_language=cpp',
    ])

add_project_arguments([ '-fno-strict-aliasing', '-w', ], language: ['c', 'cpp'])

subdir('depend-qlat')
subdir('depend-grid')

py_mod = import('python')
py3 = py_mod.find_installation('python3', pure: false)
message(py3.full_path())
message(py3.get_install_dir())

deps = [ qlat, grid, ]

subdir('qlat_grid_config')
subdir('qlat_grid')
