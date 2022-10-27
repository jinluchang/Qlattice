conf_data = configuration_data()

if cpp.has_function('mallopt', prefix: '#include <malloc.h>')
  conf_data.set('QLAT_USE_MALLOPT', true)
endif

if cpp.has_function('sysinfo', prefix: '#include <sys/sysinfo.h>')
  conf_data.set('QLAT_USE_SYSINFO', true)
endif

config = configure_file(output: 'config.h', configuration: conf_data)

install_headers(config, subdir: 'qlat')
