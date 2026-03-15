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

cmake_minimum_required(VERSION 3.17)

project(
  "gfn0"
  LANGUAGES "Fortran" "C"
  VERSION "0.1"
  DESCRIPTION "Standalone GFN0-xTB implementation"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
set(libs)
add_subdirectory("config")

###############################################################################
####################### SUBPROJECTS & DEPENDENCIES ############################
###############################################################################

# OpenMP parallelization
if(NOT TARGET "OpenMP::OpenMP_Fortran")
  if(WITH_OpenMP)
    find_package("OpenMP" REQUIRED)
    list(APPEND libs OpenMP::OpenMP_Fortran)
  else()
    add_library("OpenMP::OpenMP_Fortran" INTERFACE IMPORTED)
  endif()
endif()

# LAPACK / BLAS
if(NOT TARGET "LAPACK::LAPACK")
  find_package(LAPACK REQUIRED)
endif()
if(NOT TARGET "BLAS::BLAS")
  find_package(BLAS REQUIRED)
endif()
list(APPEND libs LAPACK::LAPACK BLAS::BLAS)

# Fortran unit test interface
if(NOT TARGET "test-drive::test-drive" AND WITH_TESTS)
  find_package("test-drive" REQUIRED)
endif()

if(WITH_GBSA)
  add_compile_definitions(WITH_GBSA)
endif()

# gfn0 sources
set(srcs)
add_subdirectory("src")

###############################################################################
############################### Static Library ################################
###############################################################################
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  WINDOWS_EXPORT_ALL_SYMBOLS TRUE
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${libs}"
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

###############################################################################
##################### Install and Export Instructions #########################
###############################################################################
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
install(
  FILES "COPYING" "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

###############################################################################
############################### Testprog / Tests ##############################
###############################################################################
if(BUILD_EXE OR WITH_TESTS)
  set(testprog-srcs)
  add_subdirectory("testprog")

  add_executable(
    "${PROJECT_NAME}-tester"
    "${testprog-srcs}"
  )
  set_target_properties(
    "${PROJECT_NAME}-tester"
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    OUTPUT_NAME "${PROJECT_NAME}-tester"
  )
  target_link_libraries(
    "${PROJECT_NAME}-tester"
    PRIVATE
    "${PROJECT_NAME}-lib"
  )
endif()

if(WITH_TESTS)
  enable_testing()
  add_test("${PROJECT_NAME}/singlepoint" "${PROJECT_NAME}-tester")
endif()
