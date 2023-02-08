init_fns = [ '__init__.py', ]
init_files = files(init_fns)

py_fns = run_command(
  'bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py | grep -v "_template\." | grep -v "__init__.py"',
  check: true,
  ).stdout().strip().split('\n')
py_files = files(py_fns)

pxd_fns = run_command(
  'bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.pxd | grep -v "_template\."',
  check: true,
  ).stdout().strip().split('\n')
pxd_files = files(pxd_fns)

pxi_fns = run_command(
  'bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.inc.pyx | grep -v "_template\."',
  check: true,
  ).stdout().strip().split('\n')
pxi_files = files(pxi_fns)

py3.install_sources(init_files + py_files + pxd_files, subdir: 'qlat_utils')

init_tgts = []
foreach i: range(init_fns.length())
  fn = init_fns[i]
  file = init_files[i]
  init_tgts += custom_target(fn,
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

pxi_tgts = []
foreach i: range(pxi_fns.length())
  fn = pxi_fns[i]
  file = pxi_files[i]
  pxi_tgts += custom_target(fn,
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
  depends : init_tgts + pxd_tgts + pxi_tgts,
  depend_files: all_header_files,
  install: false,
  )

cplib = py3.extension_module('cp',
  cp_pyx,
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat_utils',
  )

cpa_pyx = custom_target(f'cpa_pyx',
  output: 'cpa.pyx',
  input: 'cpa.pyx',
  command: [ 'cp', '@INPUT@', '@OUTPUT@', ],
  depends : init_tgts + pxd_tgts + pxi_tgts,
  depend_files: all_header_files,
  install: false,
  )

cpalib = py3.extension_module('cpa',
  cpa_pyx,
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat_utils',
  )
