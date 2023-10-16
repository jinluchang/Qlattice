conf_data = configuration_data()

if not cpp.has_function('mallopt', prefix: '#include <malloc.h>')
  conf_data.set('QLAT_NO_MALLOPT', true)
endif

if not cpp.has_function('sysinfo', prefix: '#include <sys/sysinfo.h>')
  conf_data.set('QLAT_NO_SYSINFO', true)
endif

config_file = configure_file(output: 'config.h', configuration: conf_data)

py3.install_sources(config_file, subdir: 'qlat/include/qlat')

header_files = files(
  'compressed-eigen-io.h',
  'contract-field.h',
  'contract-hvp.h',
  'contract-pion.h',
  'contract-wall-src-prop.h',
  'core.h',
  'cps.h',
  'cps-lanc.h',
  'cps-luchang.h',
  'cps-pio.h',
  'cps-utils.h',
  'dslash.h',
  'env.h',
  'fermion-action.h',
  'field-dist-io.h',
  'field-double.h',
  'field-expand.h',
  'field-fft.h',
  'field.h',
  'field-io.h',
  'field-rng.h',
  'field-serial-io.h',
  'field-shuffle.h',
  'fields-io.h',
  'flowed-hmc.h',
  'gauge-action.h',
  'geometry.h',
  'grid.h',
  'hmc.h',
  'hmc-stats.h',
  'mpi.h',
  'py_convert.h',
  'qcd-gauge-transformation.h',
  'qcd.h',
  'qcd-prop.h',
  'qcd-smear.h',
  'qcd-topology.h',
  'qcd-utils.h',
  'qed.h',
  'qlat-analysis.h',
  'qlat.h',
  'qlat-setup.h',
  'quda.h',
  'scalar-action.h',
  'selected-field.h',
  'selected-field-io.h',
  'selected-points.h',
  'utils-coordinate.h',
  'utils-io.h',
  'wilson-flow.h',
)
py3.install_sources(header_files, subdir: 'qlat/include/qlat')

py3.install_sources('qlat-setup.h', subdir: 'qlat/include')

subdir('vector_utils')

all_header_files = header_files + config_file + vector_utils_header_files
