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

cpp = meson.get_compiler('cpp')

omp = dependency('openmp')
mpic = dependency('mpi', language: 'cpp')

qlat = dependency('qlat')

deps = [ omp, mpic, qlat, ]

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/*.cpp', check: true)
sources = c.stdout().strip().split('\n')

executable('qlat.x',
  sources,
  dependencies: deps,
  install: true,
  )

run_target('run',
  command: 'run.sh')
