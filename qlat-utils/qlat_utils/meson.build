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
  'qar_utils.py',
  'qplot.py',
  'utils_io.py',
  'utils.py',
  ]
py_files = files(py_fns)

pxd_fns = [
  '__init__.pxd',
  'everything.pxd',
  'timer.pxd',
  'types.pxd',
  'mat.pxd',
  'rng_state.pxd',
  'coordinate.pxd',
  'lat_data.pxd',
  ]
pxd_files = files(pxd_fns)

py3.install_sources(init_files + py_files + pxd_files, subdir: 'qlat_utils')

pyx_dep_tgts = []
foreach f: init_files + pxd_files
  pyx_dep_tgts += fs.copyfile(f)
endforeach

py_tgts = []
foreach f: py_files
  py_tgts += fs.copyfile(f)
endforeach

pyx_name_list = [
  'timer',
  'types',
  'cutils',
  'rng_state',
  'coordinate',
  'lat_data',
  'qar',
  'mat',
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
