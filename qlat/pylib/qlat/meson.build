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

py3.install_sources(init_files + py_files + pxd_files, subdir: 'qlat')

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

type_info_list = [
  [ 'ColorMatrix', 'color_matrix', ],
  [ 'WilsonMatrix', 'wilson_matrix', ],
  [ 'NonRelWilsonMatrix', 'non_rel_wilson_matrix', ],
  [ 'IsospinMatrix', 'isospin_matrix', ],
  [ 'SpinMatrix', 'spin_matrix', ],
  [ 'WilsonVector', 'wilson_vector', ],
  [ 'Complex', 'complex', ],
  [ 'ComplexF', 'complexf', ],
  [ 'Double', 'double', ],
  [ 'Float', 'float', ],
  [ 'Long', 'long', ],
  [ 'Int64t', 'int64t', ],
  [ 'Int8t', 'int8t', ],
  [ 'Char', 'char', ],
  ]
cp_type_gen = files('cp_type_gen.sh')

cp_type_pxd_list = []
cp_type_pyx_list = []
cp_type_pxd_fn_list = []
cp_type_pyx_fn_list = []
cp_type_name_list = []

foreach ti: type_info_list
  name = ti[0]
  fn = ti[1]
  cp_type_pxd = custom_target(f'cp_type_@fn@_pxd',
    output: f'cp_type_@fn@_pxd.pxi',
    input: 'cp_type_template.pxd',
    command: [ cp_type_gen, f'@name@', '@INPUT@', '@OUTPUT@', ],
    depend_files: cp_type_gen,
    install: true,
    install_dir: py3.get_install_dir() / 'qlat',
    )
  cp_type_pyx = custom_target(f'cp_type_@fn@_pyx',
    output: f'cp_type_@fn@_pyx.pxi',
    input: 'cp_type_template.pyx',
    command: [ cp_type_gen, f'@name@', '@INPUT@', '@OUTPUT@', ],
    depend_files: cp_type_gen,
    install: false,
    )
  cp_type_pxd_list += cp_type_pxd
  cp_type_pyx_list += cp_type_pyx
  cp_type_pxd_fn_list += f'cp_type_@fn@_pxd.pxi'
  cp_type_pyx_fn_list += f'cp_type_@fn@_pyx.pxi'
  cp_type_name_list += name
endforeach

cp_gen = files('cp_gen.sh')

cp_pxd = custom_target(f'cp_pxd',
  output: 'cp.pxd',
  input: 'cp_template.pxd',
  command: [ cp_gen, '@INPUT@', '@OUTPUT@', ] + cp_type_pxd_fn_list,
  depends: cp_type_pxd_list,
  depend_files: cp_gen,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat',
  )

cp_pyx = custom_target(f'cp_pyx',
  output: 'cp.pyx',
  input: 'cp_template.pyx',
  command: [ cp_gen, '@INPUT@', '@OUTPUT@', ] + cp_type_pyx_fn_list,
  depends: cp_type_pyx_list + cp_pxd + init_tgts + pxd_tgts,
  depend_files: cp_gen + all_header_files + qlat_utils_pxd + qlat_utils_header,
  install: false,
  )

cplib = py3.extension_module('cp',
  cp_pyx,
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat',
  )

elem_type_gen = files('elem_type_gen.sh')

elem_type = custom_target(f'elem_type_py',
  output: 'elem_type.py',
  input: 'elem_type_template.py',
  command: [ elem_type_gen, '@INPUT@', '@OUTPUT@', ] + cp_type_name_list,
  depend_files: elem_type_gen,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat',
  )
