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

pxd_tgts = []
foreach i: range(pxd_fns.length())
  fn = pxd_fns[i]
  file = pxd_files[i]
  pxd_tgts += custom_target(fn,
    output: fn,
    input: file,
    command: [ 'cp', '@INPUT@', '@OUTPUT@', ],
    install: false,
    )
endforeach

py_tgts = []
foreach i: range(py_fns.length())
  fn = py_fns[i]
  file = py_files[i]
  py_tgts += custom_target(fn,
    output: fn,
    input: file,
    command: [ 'cp', '@INPUT@', '@OUTPUT@', ],
    install: false,
    )
endforeach

cp_pyx = custom_target(f'cp_pyx',
  output: 'cp.pyx',
  input: 'cp.pyx',
  command: [ 'cp', '@INPUT@', '@OUTPUT@', ],
  depends : py_tgts + pxd_tgts,
  depend_files: all_header_files + qlat_pxd + qlat_header,
  install: false,
  )

cplib = py3.extension_module('cp',
  files('cp.pyx'),
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat_cps',
  install_rpath: ':'.join(qlat_lib + [ '$ORIGIN/lib', ]),
  )
