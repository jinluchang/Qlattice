subdir('include')
subdir('lib')

py_files = files(
  '__init__.py',
  '__main__.py',
  'c.py',
  'init.py',
  'prop.py',
  )

pxd_files = files(
  '__init__.pxd',
  'everything.pxd',
  'all.pxd',
  'cp.pxd',
  )

py3.install_sources(py_files + pxd_files, subdir: 'qlat_cps')

cplib = py3.extension_module('cp',
  files('cp.pyx'),
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat_cps',
  install_rpath: ':'.join(qlat_lib + [ '$ORIGIN/lib', ]),
  )
