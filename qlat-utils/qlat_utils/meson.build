subdir('include')
subdir('lib')
subdir('scripts')

init_files = files(
  '__init__.py',
  )

py_fns = [
  'ama.py',
  'cache.py',
  'c.py',
  'data.py',
  'get_include_dir.py',
  'lat_io.py',
  'load_prop.py',
  '__main__.py',
  'parallel.py',
  'qar.py',
  'qplot.py',
  'rng_state.py',
  'utils_io.py',
  'utils.py',
  ]
py_files = files(py_fns)

pxd_fns = [
  '__init__.pxd',
  'cpa.pxd',
  'everything.pxd',
  'timer.pxd',
  'types.pxd',
  ]
pxd_files = files(pxd_fns)

pxi_fns = [
  'coordinate_d.inc.pyx',
  'coordinate.inc.pyx',
  'lat_io.inc.pyx',
  'qar.inc.pyx',
  'rng_state.inc.pyx',
  'utils.inc.pyx',
  'utils_io.inc.pyx',
  ]
pxi_files = files(pxi_fns)

py3.install_sources(init_files + py_files + pxd_files, subdir: 'qlat_utils')

pyx_dep_tgts = []
foreach f: init_files + pxd_files + pxi_files
  pyx_dep_tgts += fs.copyfile(f)
endforeach

py_tgts = []
foreach f: py_files
  py_tgts += fs.copyfile(f)
endforeach

pyx_name_list = [
  'timer',
  'types',
  'cpa',
  ]

foreach pyx_name: pyx_name_list
  pyx = custom_target(f'@pyx_name@_pyx',
    output: f'@pyx_name@.pyx',
    input: f'@pyx_name@.pyx',
    command: [ 'cp', '@INPUT@', '@OUTPUT@', ],
    depends : pyx_dep_tgts,
    depend_files: all_header_files,
    install: false,
    )
  lib = py3.extension_module(pyx_name,
    pyx,
    dependencies: deps,
    include_directories: incdir,
    link_with: clib,
    install: true,
    subdir: 'qlat_utils',
    install_rpath: '$ORIGIN/lib',
    )
endforeach
