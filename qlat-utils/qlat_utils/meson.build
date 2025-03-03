subdir('include')
subdir('lib')
subdir('scripts')

tempita = files('tempita.py')

init_files = files(
  '__init__.py',
  )

py_files = files(
  'ama.py',
  'cache.py',
  'c.py',
  'data.py',
  'load_prop.py',
  'lru_cache.py',
  '__main__.py',
  'parallel.py',
  'qplot.py',
  'utils_io.py',
  'utils.py',
  'json.py',
  'q_fit_corr.py',
  'q_fit_corr_2.py',
  'q_hlt_reconstruction.py',
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
  'qar.pxd',
  )

py3.install_sources(init_files + py_files + pxd_files, subdir: 'qlat_utils')

name_pyx_in_list = [
  'lat_data',
  ]

pxd_tgts = []
foreach f: init_files + pxd_files
  pxd_tgts += fs.copyfile(f)
endforeach

foreach name: name_pyx_in_list
  cp_pxd = custom_target(f'@name@.pxd',
    input: f'@name@.pxd.in',
    output: f'@name@.pxd',
    command: [ py3, tempita, '@INPUT@', '@OUTPUT@', ],
    install: true,
    install_dir: py3.get_install_dir() / 'qlat_utils',
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
    subdir: 'qlat_utils',
    install_rpath: '$ORIGIN/lib',
    )
endforeach

pyx_name_list = [
  'timer',
  'types',
  'cutils',
  'rng_state',
  'coordinate',
  'qar',
  'mat',
  ]

foreach pyx_name: pyx_name_list
  pyx = files(f'@pyx_name@.pyx')
  cpyx = fs.copyfile(pyx)
  lib = py3.extension_module(pyx_name,
    cpyx, pxd_tgts,
    dependencies: deps,
    include_directories: incdir,
    link_with: clib,
    install: true,
    subdir: 'qlat_utils',
    install_rpath: '$ORIGIN/lib',
    )
endforeach
