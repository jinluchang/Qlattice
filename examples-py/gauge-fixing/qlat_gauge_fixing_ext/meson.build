subdir('include')
subdir('lib')

py_files = files(
  '__init__.py',
  'c.py',
  )

pxd_files = files(
  '__init__.pxd',
  'everything.pxd',
  )

py3.install_sources(py_files + pxd_files, subdir: 'qlat_gauge_fixing_ext')

cplib = py3.extension_module('gauge_fixing',
  files('gauge_fixing.pyx'),
  dependencies: deps,
  include_directories: incdir,
  install: true,
  subdir: 'qlat_gauge_fixing_ext',
  install_rpath: ':'.join(qlat_lib + [ '$ORIGIN/lib', ]),
  )
