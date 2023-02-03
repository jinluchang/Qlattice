qlat_utils_cpp = meson.get_compiler('cpp')

qlat_utils_py3 = import('python').find_installation('python3')
message(qlat_utils_py3.path())
message(qlat_utils_py3.get_install_dir())

qlat_utils_omp = dependency('openmp').as_system()
qlat_utils_zlib = dependency('zlib').as_system()

qlat_utils_math = qlat_utils_cpp.find_library('m')

qlat_utils_numpy_include = run_command(qlat_utils_py3, '-c', 'import numpy as np ; print(np.get_include())',
  check: true).stdout().strip()
message('numpy include', qlat_utils_numpy_include)

qlat_utils_numpy = declare_dependency(
  include_directories:  include_directories(qlat_utils_numpy_include),
  dependencies: [ qlat_utils_py3.dependency(), ],
  ).as_system()

if qlat_utils_cpp.check_header('Grid/Eigen/Eigen')
  qlat_utils_eigen = dependency('', required: false)
elif qlat_utils_cpp.check_header('Eigen/Eigen')
  qlat_utils_eigen = dependency('', required: false)
else
  qlat_utils_eigen = dependency('eigen3').as_system()
endif

qlat_utils_include = run_command(qlat_utils_py3, '-c', 'import qlat_utils as q ; print("\\n".join(q.get_include_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat_utils include', qlat_utils_include)

qlat_utils_lib = run_command(qlat_utils_py3, '-c', 'import qlat_utils as q ; print("\\n".join(q.get_lib_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat_utils lib', qlat_utils_lib)

qlat_utils_pxd = run_command(qlat_utils_py3, '-c', 'import qlat_utils as q ; print("\\n".join(q.get_pxd_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat_utils pxd', qlat_utils_pxd)
qlat_utils_pxd = files(qlat_utils_pxd)

qlat_utils_header = run_command(qlat_utils_py3, '-c', 'import qlat_utils as q ; print("\\n".join(q.get_header_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat_utils header', qlat_utils_header)
qlat_utils_header = files(qlat_utils_header)

qlat_utils = declare_dependency(
  include_directories: include_directories(qlat_utils_include),
  dependencies: [
    qlat_utils_py3.dependency().as_system(),
    qlat_utils_cpp.find_library('qlat-utils', dirs: qlat_utils_lib),
    qlat_utils_numpy, qlat_utils_eigen, qlat_utils_omp, qlat_utils_zlib, qlat_utils_math, ],
  )
