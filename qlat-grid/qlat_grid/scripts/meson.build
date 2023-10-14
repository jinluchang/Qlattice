py_fns = [
  'qlat_grid_config.py',
  ]
py_files = files(py_fns)

py3.install_sources(py_files, subdir: 'qlat_grid/scripts')

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
