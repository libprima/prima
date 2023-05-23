
set(CMAKE_C_STANDARD 99)

add_library (primac cintrf.f90 cobyla_c.f90 lincoa_c.f90 bobyqa_c.f90 newuoa_c.f90 uobyqa_c.f90 prima.c)

target_include_directories (primac PUBLIC
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_link_libraries (primac PUBLIC primaf) # must be PUBLIC for precision macros
set_target_properties(primac PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 99)

install (TARGETS primac DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (FILES include/prima/prima.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prima)

macro (prima_add_c_test name)
  add_executable (example_${name}_c EXCLUDE_FROM_ALL examples/${name}/${name}_example.c)
  target_link_libraries (example_${name}_c PRIVATE primac)
  target_include_directories (example_${name}_c PRIVATE ${CMAKE_SOURCE_DIR}/c/include)
  if (PRIMA_ENABLE_EXAMPLES)
    set_target_properties (example_${name}_c PROPERTIES EXCLUDE_FROM_ALL FALSE)
  endif ()
  add_test (NAME example_${name}_c COMMAND example_${name}_c)
  add_dependencies(examples example_${name}_c)
endmacro ()

prima_add_c_test (cobyla)
prima_add_c_test (bobyqa)
prima_add_c_test (newuoa)
prima_add_c_test (uobyqa)
prima_add_c_test (lincoa)

add_subdirectory(tests)
