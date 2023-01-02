c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.py | grep -v "_template\."', check: true)
pysources = files(c.stdout().strip().split('\n'))
c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.pxd | grep -v "_template\."', check: true)
pxd_sources = files(c.stdout().strip().split('\n'))
py3.install_sources(pysources + pxd_sources, subdir: 'qlat_utils')

copy = []
foreach f: pysources
  copy += fs.copyfile(f)
endforeach

copy_pxd = []
foreach f: pxd_sources
  copy_pxd += fs.copyfile(f)
endforeach

cp_pyx = custom_target(f'cp_pyx',
  output: 'cp.pyx',
  input: 'cp.pyx',
  command: [ 'cp', '-v', '@INPUT@', '@OUTPUT@', ],
  depend_files: pxd_sources,
  install: false,
  )

cplib = py3.extension_module('cp',
  cp_pyx,
  dependencies: deps,
  include_directories: incdir,
  link_with: clib,
  install: true,
  subdir: 'qlat_utils',
  )
