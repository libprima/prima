
set(CMAKE_C_STANDARD 99)

include_directories(include)

if (WIN32 AND BUILD_SHARED_LIBS)
  add_compile_definitions(PRIMAC_EXPORTS)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    get_filename_component(Fortran_COMPILER_PATH ${CMAKE_Fortran_COMPILER} DIRECTORY)
    link_directories(${Fortran_COMPILER_PATH}/../../compiler/lib/intel64) # works for oneAPI 2023
    link_libraries(ifmodintr)
  endif ()
endif ()

get_property(PRIMAF_MODULE_DIR TARGET primaf PROPERTY Fortran_MODULE_DIRECTORY)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../fortran/${PRIMAF_MODULE_DIR}/${CMAKE_BUILD_TYPE})

add_library(primac_f OBJECT cintrf.f90 cobyla_c.f90 lincoa_c.f90 bobyqa_c.f90 newuoa_c.f90 uobyqa_c.f90)
add_library(primac_c OBJECT prima.c)
add_library(primac $<TARGET_OBJECTS:primac_f> $<TARGET_OBJECTS:primac_c>)
add_dependencies(primac_f primaf)
add_dependencies(primac_c primaf)

if (WIN32)
  set_target_properties(primac PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()

target_include_directories (primac PUBLIC
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_link_libraries (primac PUBLIC primaf) # must be PUBLIC for precision macros
set_target_properties(primac PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 99)

if (NOT BUILD_SHARED_LIBS)
  target_compile_definitions(primac PUBLIC PRIMAC_STATIC)
  target_link_libraries (primac INTERFACE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif ()

# Export symbols on Windows. See more comments in fortran/CMakeLists.txt.
if (WIN32 AND BUILD_SHARED_LIBS)
  target_sources(primac PRIVATE primac.def)
endif ()

install (FILES include/prima/prima.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prima)

macro (prima_add_c_test name)
  add_executable (example_${name}_c_exe EXCLUDE_FROM_ALL examples/${name}/${name}_example.c)
  target_link_libraries (example_${name}_c_exe PRIVATE primac)
  target_include_directories (example_${name}_c_exe PRIVATE ${CMAKE_SOURCE_DIR}/c/include)
  if (PRIMA_ENABLE_EXAMPLES)
    set_target_properties (example_${name}_c_exe PROPERTIES EXCLUDE_FROM_ALL FALSE)
  endif ()
  if (WIN32)
    set_target_properties(example_${name}_c_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
  endif()
  add_test (NAME example_${name}_c COMMAND example_${name}_c_exe)
  add_dependencies(examples example_${name}_c_exe)
endmacro ()

prima_add_c_test (cobyla)
prima_add_c_test (bobyqa)
prima_add_c_test (newuoa)
prima_add_c_test (uobyqa)
prima_add_c_test (lincoa)

add_subdirectory(tests)
