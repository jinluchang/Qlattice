project('qlat-cpp', 'cpp',
  version: '0.1',
  license: 'GPL-3.0-or-later',
  default_options: [
    'warning_level=0',
    'cpp_std=c++14',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    ])

add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])

subdir('depend-qlat')

cxx = run_command('bash', '-c', 'echo "$CXX"', check: true).stdout().strip()
mpicxx = run_command('bash', '-c', 'echo "$MPICXX"', check: true).stdout().strip()

if cxx != '' and mpicxx == cxx
  message(f'cxx=\'@cxx@\' (use CXX compiler without additional MPI options.)')
  mpic = dependency('', required: false)
else
  message(f'cxx=\'@cxx@\' mpicxx=\'@mpicxx@\' (use meson\'s automatic MPI detection.)')
  mpic = dependency('mpi', language: 'cpp').as_system()
endif

deps = [ mpic, qlat, ]

cpp_sources = run_command('bash', '-c', 'cd "$MESON_SOURCE_ROOT/$MESON_SUBDIR" ; ls *.cpp', check: true).stdout().strip().split('\n')

qlat_x = executable('qlat.x',
  cpp_sources,
  dependencies: deps,
  install: true,
  )

run_target('run',
  command: 'run.sh',
  depends: [ qlat_x, ],
  )
