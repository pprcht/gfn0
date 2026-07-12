# This file is part of gfn0.
# SPDX-Identifier: LGPL-3.0-or-later
#
# gfn0 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gfn0 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with gfn0.  If not, see <https://www.gnu.org/licenses/>.

project(
  'gfn0',
  ['fortran', 'c'],
  version         : '0.1.1',
  license         : 'LGPL-3.0-or-later',
  meson_version   : '>=0.59.0',
  default_options : [
    'buildtype=release',
    'default_library=static',
    'c_std=gnu11',
  ],
)

# ── Compiler flags and static link args ───────────────────────────────────────
# Sets: fc, cc, fc_id, cc_id, static_build
subdir('config')

inc_dirs = include_directories('include', '.')

# ═══════════════════════════════════════════════════════════════════════════════
#  OpenMP
# ═══════════════════════════════════════════════════════════════════════════════
# ifx (intel-llvm) triggers a Meson bug in dependency('openmp') because its
# preprocessor output lacks the delimiter strings Meson searches for.
# Work around it by constructing the dependency manually for that compiler.
if fc_id == 'intel-llvm'
  if get_option('openmp')
    omp_dep = declare_dependency(
      compile_args : ['-qopenmp'],
      link_args    : ['-qopenmp'],
    )
  else
    omp_dep = declare_dependency()
  endif
else
  omp_dep = dependency('openmp',
    required : get_option('openmp'),
    language : 'fortran',
  )
endif

_omp_link_dep = omp_dep.found() ? [omp_dep] : []

# ═══════════════════════════════════════════════════════════════════════════════
#  LAPACK / BLAS
# ═══════════════════════════════════════════════════════════════════════════════
lapack_opt    = get_option('lapack')
with_mkl      = false
with_openblas = false
lapack_dep    = []
blas_dep      = []

_prefer_mkl = (fc_id in ['intel', 'intel-llvm']) and (lapack_opt == 'auto')

# ── MKL ───────────────────────────────────────────────────────────────────────
if lapack_opt == 'mkl' or _prefer_mkl
  _mkl = dependency('', required : false)   # not-found sentinel

  # For static builds mkl-sdl only offers libmkl_rt.so; query the lib dir via
  # pkg-config run_command (no dependency() call) and find the real .a archives.
  if static_build
    _pc = find_program('pkg-config', required : false)
    _mkl_libdir = ''
    if _pc.found()
      _r = run_command(_pc, '--variable=libdir', 'mkl-sdl', check : false)
      if _r.returncode() == 0
        _mkl_libdir = _r.stdout().strip()
      endif
    endif
    if _mkl_libdir != ''
      _thread_lib   = (fc_id == 'gcc') ? 'mkl_gnu_thread' : 'mkl_intel_thread'
      _mkl_sa_parts = []
      foreach lib : ['mkl_intel_lp64', _thread_lib, 'mkl_core']
        _l = fc.find_library(lib, dirs : [_mkl_libdir], required : false, static : true)
        if _l.found()
          _mkl_sa_parts += [_l]
        endif
      endforeach
      if _mkl_sa_parts.length() == 3
        _mkl = declare_dependency(
          dependencies : _mkl_sa_parts,
          link_args    : ['-lpthread', '-lm', '-ldl'],
        )
      endif
    endif
  endif

  if not _mkl.found()
    _mkl = dependency('mkl-sdl', required : false, static : static_build)
  endif

  if not _mkl.found()
    _thread_lib = (fc_id == 'gcc') ? 'mkl_gnu_thread' : 'mkl_intel_thread'
    _omp_lib    = (fc_id == 'gcc') ? 'gomp'           : 'iomp5'
    _mkl_parts  = []
    foreach lib : ['mkl_intel_lp64', _thread_lib, 'mkl_core', _omp_lib, 'm', 'dl']
      _l = fc.find_library(lib, required : false, static : static_build)
      if _l.found()
        _mkl_parts += [_l]
      endif
    endforeach
    if _mkl_parts.length() > 0
      _mkl = declare_dependency(dependencies : _mkl_parts)
    endif
  endif
  if _mkl.found()
    lapack_dep = _mkl
    blas_dep   = _mkl
    with_mkl   = true
    message('LAPACK/BLAS provider: Intel MKL')
  elif lapack_opt == 'mkl'
    error('MKL requested but not found.')
  else
    lapack_opt = 'auto'
  endif
endif

# ── OpenBLAS ──────────────────────────────────────────────────────────────────
if lapack_opt in ['auto', 'openblas'] and not with_mkl
  _opa = dependency('openblas', required : false, static : static_build)
  if not _opa.found()
    _opa = dependency('blas', required : false, static : static_build)
  endif
  if _opa.found()
    _opa_has_lapack = fc.links(
      'external dsygvd; call dsygvd(); end',
      dependencies : _opa,
    )
    lapack_dep    = _opa_has_lapack ? _opa : dependency('lapack',
                      required : false, static : static_build)
    blas_dep      = _opa
    with_openblas = true
    message('LAPACK/BLAS provider: OpenBLAS')
  elif lapack_opt == 'openblas'
    error('OpenBLAS requested but not found.')
  endif
endif

# ── Netlib reference ──────────────────────────────────────────────────────────
if lapack_opt in ['auto', 'netlib'] and not with_mkl and not with_openblas
  _blas   = dependency('blas',   required : false, static : static_build)
  _lapack = dependency('lapack', required : false, static : static_build)
  if _blas.found() and _lapack.found()
    blas_dep   = _blas
    lapack_dep = _lapack
    message('LAPACK/BLAS provider: netlib reference')
  else
    _blas   = fc.find_library('blas',   required : false, static : static_build)
    _lapack = fc.find_library('lapack', required : false, static : static_build)
    if _blas.found() and _lapack.found()
      blas_dep   = declare_dependency(dependencies : [_blas])
      lapack_dep = declare_dependency(dependencies : [_lapack])
    elif lapack_opt == 'netlib'
      error('Netlib BLAS/LAPACK not found.')
    else
      error('No LAPACK/BLAS found. Install libopenblas-dev or set -Dlapack=custom.')
    endif
  endif
endif

# ── Custom ────────────────────────────────────────────────────────────────────
if lapack_opt == 'custom'
  _ll = get_option('lapack_libs')
  _bl = get_option('blas_libs')
  if _ll.length() == 0 and _bl.length() == 0
    error('lapack=custom requires lapack_libs and/or blas_libs.')
  endif
  _deps = []
  foreach lib : _ll + _bl
    _deps += [fc.find_library(lib, required : true, static : static_build)]
  endforeach
  lapack_dep = declare_dependency(dependencies : _deps)
  blas_dep   = lapack_dep
endif

if with_mkl
  add_project_arguments('-DWITH_MKL',      language : ['c', 'fortran'])
endif
if with_openblas
  add_project_arguments('-DWITH_OPENBLAS', language : ['c', 'fortran'])
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Optional features
# ═══════════════════════════════════════════════════════════════════════════════
with_gbsa = get_option('gbsa').allowed()
if with_gbsa
  add_project_arguments('-DWITH_GBSA', language : 'fortran')
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Source collection
# ═══════════════════════════════════════════════════════════════════════════════
srcs = []
subdir('src')

# ═══════════════════════════════════════════════════════════════════════════════
#  Library target
# ═══════════════════════════════════════════════════════════════════════════════
_lib_deps = [lapack_dep, blas_dep] + _omp_link_dep

gfn0_lib = library(
  meson.project_name(),
  sources             : srcs,
  include_directories : inc_dirs,
  dependencies        : _lib_deps,
  install             : not (meson.is_subproject() and get_option('default_library') == 'static'),
)

# ── Exported dependency (what crest links against) ────────────────────────────
gfn0_dep = declare_dependency(
  link_with           : gfn0_lib,
  include_directories : [inc_dirs, gfn0_lib.private_dir_include()],
  dependencies        : _lib_deps,
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Testprog — standalone executable, optionally registered as a test
# ═══════════════════════════════════════════════════════════════════════════════
# The testprog is not a test-drive suite but a self-contained program that
# runs several GFN0 API calls on a hardcoded caffeine molecule and exits 0
# on success.  build_exe=true builds it; tests=true additionally registers
# it with `meson test`.  You can also set build_exe=false, tests=true —
# the executable is still built in that case (it is needed to run the test).
if get_option('build_exe') or get_option('tests')
  subdir('testprog')
  gfn0_tester = executable(
    'gfn0-tester',
    sources             : prog,
    include_directories : inc_dirs,
    dependencies        : gfn0_dep,
    link_language       : 'fortran',
    install             : get_option('build_exe'),
  )
  if get_option('tests')
    test('gfn0/singlepoint', gfn0_tester, timeout : 120)
  endif
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Install
# ═══════════════════════════════════════════════════════════════════════════════
install_data(
  'LICENSE',
  install_dir : get_option('datadir') / 'licenses' / meson.project_name(),
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════════════════════
summary({
  'Fortran compiler' : '@0@ @1@'.format(fc_id, fc.version()),
  'C compiler'       : '@0@ @1@'.format(cc_id, cc.version()),
  'Build type'       : get_option('buildtype'),
  'Static binary'    : static_build,
  'OpenMP'           : omp_dep.found(),
  'MKL'              : with_mkl,
  'OpenBLAS'         : with_openblas,
  'GBSA/ALPB'        : with_gbsa,
  'Testprog'         : get_option('build_exe'),
  'Tests registered' : get_option('tests'),
}, section : 'gfn0 build configuration')
