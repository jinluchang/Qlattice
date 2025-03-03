py_fns = [
  'eigen_system_checksum.py',
  'eigen_system_repartition.py',
  'fields_checksum.py',
  'fields_list.py',
  'fields_build_index.py',
  'fields_properly_truncate.py',
  'fields_rewrite.py',
  'gauge_fix_coulomb.py',
  'topo_measure.py',
  'topo_measure_gpt.py',
  'topo_measure_wilson_flow.py',
  ]
py_files = files(py_fns)

py3.install_sources(py_files, subdir: 'qlat/scripts')

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
