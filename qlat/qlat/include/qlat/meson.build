conf_data = configuration_data()

if not cpp.has_function('mallopt', prefix: '#include <malloc.h>')
  conf_data.set('QLAT_NO_MALLOPT', true)
endif

if not cpp.has_function('sysinfo', prefix: '#include <sys/sysinfo.h>')
  conf_data.set('QLAT_NO_SYSINFO', true)
endif

if not cpp.has_function('Cuhre', prefix: '#include <cuba.h>', dependencies: cuba)
  conf_data.set('QLAT_NO_CUBA', true)
endif

config_file = configure_file(output: 'config.h', configuration: conf_data)

py3.install_sources(config_file, subdir: 'qlat/include/qlat')

header_files = files(
  'compressed-eigen-io.h',
  'compute-f.h',
  'compute-int-mult.h',
  'compute-int-seq.h',
  'contract-field.h',
  'contract-hvp.h',
  'contract-pion.h',
  'contract-wall-src-prop.h',
  'core.h',
  'cps-lanc.h',
  'cps-luchang.h',
  'cps-pio.h',
  'cps-utils.h',
  'cps.h',
  'dslash.h',
  'env.h',
  'fermion-action.h',
  'field-base-io.h',
  'field-dist-io.h',
  'field-double.h',
  'field-expand.h',
  'field-fft.h',
  'field-io.h',
  'field-rng.h',
  'field-serial-io.h',
  'field-shuffle.h',
  'field.h',
  'fields-io.h',
  'flowed-hmc.h',
  'function.h',
  'gauge-action.h',
  'geometry.h',
  'grid.h',
  'grid-sys.h',
  'hmc-stats.h',
  'hmc.h',
  'integration-multidimensional.h',
  'integration.h',
  'interpolation-bilinear.h',
  'interpolation.h',
  'mpi.h',
  'muon-line-config.h',
  'muon-line.h',
  'hlbl-sl-table.h',
  'hlbl-contract.h',
  'projection.h',
  'py_convert.h',
  'qcd-gauge-transformation.h',
  'qcd-prop.h',
  'qcd-smear.h',
  'qcd-topology.h',
  'qcd-utils.h',
  'qcd-acc.h',
  'qcd.h',
  'qed.h',
  'qlat-analysis.h',
  'qlat-setup.h',
  'qlat.h',
  'qm-action.h',
  'quda.h',
  'root-fsolver.h',
  'scalar-action.h',
  'selected-field-io.h',
  'selected-field.h',
  'selected-points.h',
  'selected-shuffle.h',
  'utils-coordinate.h',
  'utils-io.h',
  'utils.h',
  'wilson-flow.h',
)
py3.install_sources(header_files, subdir: 'qlat/include/qlat')

py3.install_sources('qlat-setup.h', subdir: 'qlat/include')

subdir('vector_utils')

all_header_files = header_files + config_file + vector_utils_header_files
