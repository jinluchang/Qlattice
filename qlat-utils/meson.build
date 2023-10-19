project('qlat-utils', 'cpp', 'cython',
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

cpp = meson.get_compiler('cpp')

py_mod = import('python')
py3 = py_mod.find_installation('python3', pure: false)
message(py3.path())
message(py3.get_install_dir())

fs = import('fs')

omp = dependency('openmp').as_system()
zlib = dependency('zlib').as_system()
math = cpp.find_library('m')

numpy_include = run_command(py3.path(), '-c', 'import numpy as np ; print(np.get_include())', check: true).stdout().strip()
message('numpy include', numpy_include)

numpy = declare_dependency(
  include_directories:  include_directories(numpy_include),
  dependencies: [ py3.dependency(), ],
  ).as_system()

if cpp.check_header('Eigen/Eigen')
  eigen = dependency('', required: false)
elif cpp.check_header('Grid/Eigen/Eigen')
  eigen = dependency('', required: false)
else
  eigen = dependency('eigen3').as_system()
endif

deps = [ py3.dependency().as_system(), numpy, eigen, omp, zlib, math, ]

subdir('qlat_utils')
subdir('bin')
