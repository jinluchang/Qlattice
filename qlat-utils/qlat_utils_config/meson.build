config = configure_file(
  capture: true,
  command: [
    'echo', f'eigen_type = "@eigen_type@"'
    ],
  output: 'config.py',
  )

init_files = files(
  '__init__.py',
  )

py_files = files(
  'get_info.py',
  '__main__.py',
  ) + config

py3.install_sources(init_files + py_files, subdir: 'qlat_utils_config')
