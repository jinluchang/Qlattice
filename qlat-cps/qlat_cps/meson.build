subdir('scripts')
subdir('include')
subdir('lib')

py_fns = [
  'c.py',
  'get_include_dir.py',
  '__init__.py',
  'init.py',
  '__main__.py',
  'prop.py',
  ]
py_files = files(py_fns)

pxd_fns = [
  '__init__.pxd',
  'cp.pxd',
  'everything.pxd',
  ]
pxd_files = files(pxd_fns)

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
