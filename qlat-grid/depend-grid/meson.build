grid_cpp = meson.get_compiler('cpp')

message('Collecting Grid information using \'grid-config\':')
grid_prefix = run_command('grid-config', '--prefix', check: true).stdout().strip()
message('Grid: prefix =', grid_prefix)
message('Grid: CXX =', run_command('grid-config', '--cxx', check: true).stdout().strip())
message('Grid: CXXLD =', run_command('grid-config', '--cxxld', check: true).stdout().strip())
grid_cxxflags = run_command('grid-config', '--cxxflags', check: true).stdout().strip().split()
grid_ldflags = run_command('grid-config', '--ldflags', check: true).stdout().strip().split()
grid_libs = run_command('grid-config', '--libs', check: true).stdout().strip().split()
message('Grid: CXXFLAGS =', grid_cxxflags)
message('Grid: LDFLAGS =', grid_ldflags)
message('Grid: LIBS =', grid_libs)

grid_compile_args = []
grid_include_directories = []
foreach f : grid_cxxflags
  if f.startswith('-I')
    v = f.substring(2)
    if v not in grid_include_directories
      grid_include_directories += v
    endif
  else
    if f not in grid_compile_args
      grid_compile_args += f
    endif
  endif
endforeach

grid_link_args = []
grid_lib_dirs = []
grid_lib_names = [ 'Grid', ]
foreach f : grid_ldflags + grid_libs
  if f.startswith('-L')
    v = f.substring(2)
    if v not in grid_lib_dirs
      grid_lib_dirs += v
    endif
  elif f.startswith('-l')
    v = f.substring(2)
    if v not in grid_lib_names
      grid_lib_names += v
    endif
  else
    grid_link_args += f
  endif
endforeach

message('Grid processed info:')
message('Grid: compile_args =', grid_compile_args)
message('Grid: include_directories =', grid_include_directories)
message('Grid: link_args =', grid_link_args)
message('Grid: lib_dirs =', grid_lib_dirs)
message('Grid: lib_names =', grid_lib_names)

grid_deps = []
grid_link_args_not_found_libs = []
foreach lib : grid_lib_names
  dep = grid_cpp.find_library(lib, dirs: grid_lib_dirs, required: false)
  if dep.found()
    grid_deps += dep
  else
    grid_link_args_not_found_libs += '-l' + lib
  endif
endforeach

message('Grid: link_args_not_found_libs =', grid_link_args_not_found_libs)

grid = declare_dependency(
  include_directories: grid_include_directories,
  compile_args: grid_compile_args,
  link_args: grid_link_args + grid_link_args_not_found_libs,
  dependencies: grid_deps,
  ).as_system()
