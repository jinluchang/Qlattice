subdir('include')

py_fns = [
  '__init__.py',
  'c.py',
  ]
py_files = files(py_fns)

pxd_fns = [
  'hello.pxd',
  'everything.pxd',
  'all.pxd',
  ]
pxd_files = files(pxd_fns)

py3.install_sources(py_files + pxd_files, subdir: 'qlat_ext')

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
  output: 'hello.pyx',
  input: 'hello.pyx',
  command: [ 'cp', '@INPUT@', '@OUTPUT@', ],
  depends : py_tgts + pxd_tgts,
  depend_files: all_header_files + qlat_pxd + qlat_header,
  install: false,
  )

cplib = py3.extension_module('hello',
  cp_pyx,
  dependencies: deps,
  include_directories: incdir,
  install: true,
  subdir: 'qlat_ext',
  # gnu_symbol_visibility: 'default',
  install_rpath: ':'.join(qlat_lib),
  )
