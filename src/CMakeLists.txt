
IF(WITH_GBSA)
add_subdirectory("gbsa")
ENDIF()

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/gfn0_math_wrapper.f90"
  "${dir}/slater_module.f90"
  "${dir}/wfn_module.f90"
  "${dir}/cn_module.f90"
  "${dir}/eeq_module.F90"
  "${dir}/dftd4param.f90"
  "${dir}/overlap_module.f90"
  "${dir}/gfn0_types.f90"
  "${dir}/gfn0_dftd4.f90"
  "${dir}/gfn0_srb.f90"
  "${dir}/gfn0_param.f90"
  "${dir}/gfn0_basis.f90"
  "${dir}/gfn0_qm.f90"
  "${dir}/gfn0_prints.f90"
  "${dir}/gfn0_gbsa.F90" 
  "${dir}/gfn0_occ.f90"
  "${dir}/gfn0_module.f90"
  "${dir}/gfn0_interface.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
