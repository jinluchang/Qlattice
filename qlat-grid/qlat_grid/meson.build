subdir('scripts')
subdir('include')
subdir('lib')

py_files = files(
  'c.py',
  'get_include_dir.py',
  '__init__.py',
  'init.py',
  '__main__.py',
  'prop.py',
  )

pxd_files = files(
  '__init__.pxd',
  'all.pxd',
  'cp.pxd',
  'everything.pxd',
  )

py3.install_sources(py_files + pxd_files, subdir: 'qlat_grid')

cplib = py3.extension_module('cp',
  files('cp.pyx'),
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat_grid',
  install_rpath: ':'.join(qlat_lib + [ '$ORIGIN/lib', ]),
  )
