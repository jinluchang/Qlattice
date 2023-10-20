subdir('scripts')
subdir('include')
subdir('lib')

tempita = files('tempita.py')

init_files = files(
  '__init__.py',
  )

py_files = files(
  'contract_field.py',
  'contract_hvp.py',
  'contract_pion.py',
  'coordinate.py',
  'fermion_action.py',
  'field_double.py',
  'field_selection_utils.py',
  'fields_io.py',
  'field_utils.py',
  'field_type_dict.py',
  'fthmc.py',
  'gauge_action.py',
  'get_include_dir.py',
  'hmc.py',
  'hmc_stats.py',
  'inverter.py',
  '__main__.py',
  'mpi.py',
  'propagator.py',
  'qcd.py',
  'scalar_action.py',
  'smear.py',
  'topology.py',
  'wilson_flow.py',
  )

pxd_files = files(
  '__init__.pxd',
  'all.pxd',
  'everything.pxd',
  'geometry.pxd',
  'field_selection.pxd',
  'field_base.pxd',
  )

py3.install_sources(init_files + py_files + pxd_files, subdir: 'qlat')

pxd_tgts = []
foreach f: init_files + pxd_files
  pxd_tgts += fs.copyfile(f)
endforeach

c_py = custom_target('c.py',
  input: 'c.py.in',
  output: 'c.py',
  command: [ py3, tempita, '@INPUT@', '@OUTPUT@', ],
  install: true,
  install_dir: py3.get_install_dir() / 'qlat',
  )

name_pyx_in_list = [
  'field_types',
  'selected_field_types',
  'selected_points_types',
  ]

foreach name: name_pyx_in_list
  cp_pxd = custom_target(f'@name@.pxd',
    input: f'@name@.pxd.in',
    output: f'@name@.pxd',
    command: [ py3, tempita, '@INPUT@', '@OUTPUT@', ],
    install: true,
    install_dir: py3.get_install_dir() / 'qlat',
    )
  pxd_tgts += cp_pxd
endforeach

foreach name: name_pyx_in_list
  cp_pyx = custom_target(f'@name@.pyx',
    input: f'@name@.pyx.in',
    output: f'@name@.pyx',
    command: [ py3, tempita, '@INPUT@', '@OUTPUT@', ],
    install: false,
    )
  cplib = py3.extension_module(f'@name@',
    cp_pyx, pxd_tgts,
    dependencies: deps,
    include_directories: incdir,
    link_with: clib,
    install: true,
    subdir: 'qlat',
    install_rpath: ':'.join(qlat_utils_lib + [ '$ORIGIN/lib', ]),
    )
endforeach

name_pyx_list = [
  'geometry',
  'field_selection',
  'field_base',
  'utils_io',
  ]

foreach name: name_pyx_list
  pyx = files(f'@name@.pyx')
  cpyx = fs.copyfile(pyx)
  lib = py3.extension_module(name,
    cpyx, pxd_tgts,
    dependencies: deps,
    include_directories: incdir,
    link_with: clib,
    install: true,
    subdir: 'qlat',
    install_rpath: ':'.join(qlat_utils_lib + [ '$ORIGIN/lib', ]),
    )
endforeach
