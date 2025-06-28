subdir('include')
subdir('lib')

py_files = files(
  '__init__.py',
  'c.py',
  )

pxd_files = files(
  '__init__.pxd',
  'everything.pxd',
  'all.pxd',
  'cp.pxd',
  )

py3.install_sources(py_files + pxd_files, subdir: 'qlat_ext')

cplib = py3.extension_module('cp',
  files('cp.pyx'),
  dependencies: deps,
  include_directories: incdir,
  install: true,
  subdir: 'qlat_ext',
  install_rpath: ':'.join(qlat_lib + [ '$ORIGIN/lib', ]),
  )
