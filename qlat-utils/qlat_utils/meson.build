subdir('include')
subdir('lib')
subdir('scripts')

py_files = files(
  '__init__.py',
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
  )

pxd_files = files(
  '__init__.pxd',
  'all.pxd',
  'everything.pxd',
  'timer.pxd',
  'types.pxd',
  'mat.pxd',
  'rng_state.pxd',
  'coordinate.pxd',
  'lat_data.pxd',
  )

py3.install_sources(py_files + pxd_files, subdir: 'qlat_utils')

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
  lib = py3.extension_module(pyx_name,
    f'@pyx_name@.pyx',
    dependencies: deps,
    include_directories: incdir,
    link_with: clib,
    install: true,
    subdir: 'qlat_utils',
    install_rpath: '$ORIGIN/lib',
    )
endforeach
