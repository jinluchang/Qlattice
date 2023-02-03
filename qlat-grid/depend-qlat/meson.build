qlat_cpp = meson.get_compiler('cpp')

qlat_py3 = import('python').find_installation('python3')
message(qlat_py3.path())
message(qlat_py3.get_install_dir())

qlat_omp = dependency('openmp').as_system()
qlat_zlib = dependency('zlib').as_system()

qlat_fftw = dependency('fftw3').as_system()
qlat_fftwf = dependency('fftw3f').as_system()
message('fftw libdir', qlat_fftw.get_variable('libdir'))
message('fftwf libdir', qlat_fftwf.get_variable('libdir'))
qlat_fftw_all = [ qlat_fftw, qlat_fftwf, ]

qlat_math = qlat_cpp.find_library('m')

qlat_numpy_include = run_command(qlat_py3, '-c', 'import numpy as np ; print(np.get_include())',
  check: true).stdout().strip()
message('numpy include', qlat_numpy_include)

qlat_numpy = declare_dependency(
  include_directories:  include_directories(qlat_numpy_include),
  dependencies: [ qlat_py3.dependency(), ],
  ).as_system()

if qlat_cpp.check_header('Grid/Eigen/Eigen')
  qlat_eigen = dependency('', required: false)
elif qlat_cpp.check_header('Eigen/Eigen')
  qlat_eigen = dependency('', required: false)
else
  qlat_eigen = dependency('eigen3').as_system()
endif

qlat_include = run_command(qlat_py3, '-c', 'import qlat as q ; print("\\n".join(q.get_include_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat include', qlat_include)

qlat_lib = run_command(qlat_py3, '-c', 'import qlat as q ; print("\\n".join(q.get_lib_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat lib', qlat_lib)

qlat_pxd = run_command(qlat_py3, '-c', 'import qlat as q ; print("\\n".join(q.get_pxd_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat pxd', qlat_pxd)
qlat_pxd = files(qlat_pxd)

qlat_header = run_command(qlat_py3, '-c', 'import qlat as q ; print("\\n".join(q.get_header_list()))',
  env: environment({'q_verbose': '-1'}),
  check: true).stdout().strip().split('\n')
message('qlat header', qlat_header)
qlat_header = files(qlat_header)

qlat = declare_dependency(
  include_directories: include_directories(qlat_include),
  dependencies: [
    qlat_py3.dependency().as_system(),
    qlat_cpp.find_library('qlat', dirs: qlat_lib),
    qlat_cpp.find_library('qlat-utils', dirs: qlat_lib),
    qlat_numpy, qlat_eigen, qlat_omp, qlat_fftw_all, qlat_zlib, qlat_math, ],
  )
