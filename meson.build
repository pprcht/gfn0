project(
  'gfn0',
  'fortran', 'c',
  version: '3.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.63',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=static',
    'c_link_args=-static',
    'fortran_link_args=-static',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
exe_deps = []
subdir('config')

# Documentation
#subdir('docs')

# Collect source of the project
prog = []
srcs = []
subdir('src')


# Create library target
gfn0_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: exe_deps,
  include_directories: include_directories('include'),
)

# Export as dependency
gfn0_inc = gfn0_lib.private_dir_include()
gfn0_dep = declare_dependency(
  link_with: gfn0_lib,
  include_directories: gfn0_inc,
  dependencies: exe_deps,
)


# Create executable target
if get_option('build_exe')
gfn0_exe = executable(
  meson.project_name(),
  sources: prog,
  dependencies: gfn0_dep,
  install: install,
  link_language: 'fortran',
)
endif


# add the testsuite separate meson.build 
#subdir('testsuite') #has to be filled with tests, not availabel yet
