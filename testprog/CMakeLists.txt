# This file is part of gfn0.
# SPDX-Identifier: LGPL-3.0-or-later

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND testprog-srcs
  "${dir}/testmol.f90"
  "${dir}/main.f90"
)

set(testprog-srcs ${testprog-srcs} PARENT_SCOPE)
