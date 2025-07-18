project('qlat-cps', 'cpp', 'cython',
  version: '0.80',
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

py_mod = import('python')
py3 = py_mod.find_installation('python3', pure: false)
message(py3.full_path())
message(py3.get_install_dir())

cxx = run_command('bash', '-c', 'echo "$CXX"', check: true).stdout().strip()
mpicxx = run_command('bash', '-c', 'echo "$MPICXX"', check: true).stdout().strip()

if cxx != '' and mpicxx == cxx
  message(f'cxx=\'@cxx@\' (use CXX compiler without additional MPI options.)')
  mpic = dependency('', required: false)
else
  message(f'cxx=\'@cxx@\' mpicxx=\'@mpicxx@\' (use meson\'s automatic MPI detection.)')
  mpic = dependency('mpi', language: 'cpp').as_system()
endif

subdir('depend-qlat')
subdir('depend-cps')

deps = [ mpic, qlat, cps, ]

subdir('qlat_cps_config')
subdir('qlat_cps')
