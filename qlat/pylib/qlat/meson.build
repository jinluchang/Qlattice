c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py', check: true)
pysources = files(c.stdout().strip().split('\n'))
c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.pxd | grep -v "cp_type_template.pxd\|cp.pxd"', check: true)
pxd_sources = files(c.stdout().strip().split('\n'))
py3.install_sources(pysources + pxd_sources, subdir: 'qlat')

copy = []
foreach f: pysources
  copy += fs.copyfile(f)
endforeach

copy_pxd = []
foreach f: pxd_sources
  copy_pxd += fs.copyfile(f)
endforeach

type_info_list = [
  [ 'ColorMatrix', 'color_matrix', ],
  [ 'WilsonMatrix', 'wilson_matrix', ],
  [ 'NonRelWilsonMatrix', 'non_rel_wilson_matrix', ],
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

cp_gen = files('cp_gen.sh')

cp_type_pxd_list = []
cp_type_pyx_list = []

cp_type_pxd_fn_list = []
cp_type_pyx_fn_list = []

foreach ti: type_info_list
  name = ti[0]
  fn = ti[1]
  cp_type_pxd = custom_target(f'cp_type_@fn@_pxd',
    output: f'cp_type_@fn@_pxd.pxi',
    input: 'cp_type_template.pxd',
    command: [ cp_type_gen, f'@name@', '@INPUT@', '@OUTPUT@', ],
    install: true,
    install_dir: py3.get_install_dir() / 'qlat',
    )
  cp_type_pyx = custom_target(f'cp_type_@fn@_pyx',
    output: f'cp_type_@fn@_pyx.pxi',
    input: 'cp_type_template.pyx',
    command: [ cp_type_gen, f'@name@', '@INPUT@', '@OUTPUT@', ],
    install: false,
    )
  cp_type_pxd_list += cp_type_pxd
  cp_type_pyx_list += cp_type_pyx
  cp_type_pxd_fn_list += f'cp_type_@fn@_pxd.pxi'
  cp_type_pyx_fn_list += f'cp_type_@fn@_pyx.pxi'
endforeach

cp_pxd = custom_target(f'cp_pxd',
  output: 'cp.pxd',
  input: 'cp.pxd',
  command: [ cp_gen, '@INPUT@', '@OUTPUT@', ] + cp_type_pxd_fn_list,
  depends: cp_type_pxd_list,
  install: true,
  install_dir: py3.get_install_dir() / 'qlat',
  )

cp_pyx = custom_target(f'cp_pyx',
  output: 'cp.pyx',
  input: 'cp.pyx',
  command: [ cp_gen, '@INPUT@', '@OUTPUT@', ] + cp_type_pyx_fn_list,
  depends: cp_type_pyx_list + cp_pxd + copy_pxd,
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
