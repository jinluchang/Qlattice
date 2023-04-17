cps_cpp = meson.get_compiler('cpp')

message('Make sure the libraries (libcps.a etc) are in the paths listed in the LIBRARY_PATH environment variable, and header files are in the paths listed in the C_INCLUDE_PATH and CPLUS_INCLUDE_PATH environment variable.')

add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])

cps_deps = []

cps_omp = dependency('openmp').as_system()
cps_fftw = dependency('fftw3').as_system()
cps_fftwf = dependency('fftw3f').as_system()
message('fftw libdir', cps_fftw.get_variable('libdir'))
message('fftwf libdir', cps_fftwf.get_variable('libdir'))
cps_fftw_all = [ cps_fftw, cps_fftwf, ]

cps_deps += [ cps_omp, cps_fftw_all, ]

cps_lib_names = [ 'cps', 'qio', 'qmp', 'lime', ]

cps_compile_args = [ '-DUSE_OMP', '-DUSE_C11', '-DUSE_QMP', '-DUSE_QIO', ]
cps_link_args_not_found_libs = []

foreach lib : cps_lib_names
  dep = cps_cpp.find_library(lib, required: false)
  if dep.found()
    cps_deps += dep
  else
    cps_link_args_not_found_libs += '-l' + lib
  endif
endforeach

message('CPS processed info:')
message('CPS: lib_names =', cps_lib_names)
message('CPS: compile_args =', cps_compile_args)
message('CPS: link_args_not_found_libs =', cps_link_args_not_found_libs)

cps = declare_dependency(
  include_directories: [],
  compile_args: cps_compile_args,
  link_args: cps_link_args_not_found_libs,
  dependencies: cps_deps,
  ).as_system()
