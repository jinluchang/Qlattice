py_fns = run_command(
  'bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py | grep -v "_template\."',
  check: true,
  ).stdout().strip().split('\n')
py_files = files(py_fns)

pxd_fns = run_command(
  'bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.pxd | grep -v "_template\."',
  check: true,
  ).stdout().strip().split('\n')
pxd_files = files(pxd_fns)

py3.install_sources(py_files + pxd_files, subdir: 'qlat_grid')

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
  install: true,
  subdir: 'qlat_grid',
  # gnu_symbol_visibility: 'default',
  )
