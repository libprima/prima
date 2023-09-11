
add_library (primaf 
  common/old/ieee_4dev.f90
  common/linalg.f90
  common/pintrf.f90
  common/string.f90
  common/evaluate.f90
  common/preproc.f90
  common/univar.f90
  common/powalg.f90
  common/history.f90
  common/xinbd.f90
  common/message.f90
  common/fprint.f90
  common/ratio.f90
  common/debug.F90
  common/consts.F90
  common/checkexit.f90
  common/inf.F90
  common/redrho.f90
  common/memory.F90
  common/selectx.f90
  common/infnan.F90
  common/shiftbase.f90
  common/infos.f90
  cobyla/cobyla.f90
  cobyla/cobylb.f90
  cobyla/geometry.f90
  cobyla/initialize.f90
  cobyla/trustregion.f90
  cobyla/update.f90
  bobyqa/bobyqa.f90
  bobyqa/bobyqb.f90
  bobyqa/geometry.f90
  bobyqa/initialize.f90
  bobyqa/trustregion.f90
  bobyqa/rescue.f90
  bobyqa/update.f90
  lincoa/update.f90
  lincoa/initialize.f90
  lincoa/getact.f90
  lincoa/trustregion.f90
  lincoa/geometry.f90
  lincoa/lincob.f90
  lincoa/lincoa.f90
  newuoa/initialize.f90
  newuoa/trustregion.f90
  newuoa/geometry.f90
  newuoa/update.f90
  newuoa/newuob.f90
  newuoa/newuoa.f90
  uobyqa/initialize.f90
  uobyqa/update.f90
  uobyqa/geometry.f90
  uobyqa/trustregion.f90
  uobyqa/uobyqb.f90
  uobyqa/uobyqa.f90
)

set (PRIMA_REAL_PRECISION "64" CACHE STRING "real precision")
set (PRIMA_INTEGER_KIND "0" CACHE STRING "integer kind")
target_compile_definitions (primaf PUBLIC "PRIMA_REAL_PRECISION=${PRIMA_REAL_PRECISION}")
target_compile_definitions (primaf PUBLIC "PRIMA_INTEGER_KIND=${PRIMA_INTEGER_KIND}")
set_target_properties(primaf PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY mod)
target_include_directories (primaf PUBLIC
  $<INSTALL_INTERFACE:include/prima/mod>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

# because of calcfc_internal nested subroutine in cobylb.f90
if (HAVE_WARN_EXECSTACK)
  target_link_options (primaf PUBLIC "-Wl,--no-warn-execstack")
endif ()

install (TARGETS primaf DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prima)

macro (prima_add_f_test name)
  add_executable (example_${name}_fortran EXCLUDE_FROM_ALL examples/${name}/${name}_example.f90)
  target_link_libraries (example_${name}_fortran PRIVATE primaf)
  target_include_directories (example_${name}_fortran PRIVATE ${CMAKE_BINARY_DIR}/fortran)
  set_target_properties(example_${name}_fortran PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/examples/${name}/mod)
  if (PRIMA_ENABLE_EXAMPLES)
    set_target_properties (example_${name}_fortran PROPERTIES EXCLUDE_FROM_ALL FALSE)
  endif ()
  add_test (NAME example_${name}_f COMMAND example_${name}_fortran)
  add_dependencies(examples example_${name}_fortran examples)
endmacro ()

prima_add_f_test (cobyla)
prima_add_f_test (bobyqa)
prima_add_f_test (newuoa)
prima_add_f_test (uobyqa)
prima_add_f_test (lincoa)
