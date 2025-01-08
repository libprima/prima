# For Apple and Windows we need to make sure it doesn't try to find the system python.
# Details: https://github.com/actions/setup-python/issues/121#issuecomment-777748504
if ((APPLE OR WIN32) AND DEFINED ENV{CI})
  if (POLICY CMP0094)  # https://cmake.org/cmake/help/latest/policy/CMP0094.html
      cmake_policy(SET CMP0094 NEW)  # FindPython should return the first matching Python
  endif ()
  set(Python_FIND_FRAMEWORK "NEVER")
  set(Python_FIND_REGISTRY "NEVER")
endif ()

set(PYBIND11_NEWPYTHON ON)
add_subdirectory(pybind11)


# numpy is needed since we use pybind11::array and array_t
find_package (Python 3.6 REQUIRED COMPONENTS Interpreter)
execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy" RESULT_VARIABLE _IMPORT_NUMPY)
if (NOT _IMPORT_NUMPY EQUAL 0)
  message(SEND_ERROR "numpy: NOT FOUND, not installing/testing Python bindings")
  return()
endif ()


# We add NO_EXTRAS to disable link time optimization. This is a workaround for musl
# due to some issue about 'inlining failed in call to always_inline vsnprintf'.
# https://gitlab.alpinelinux.org/alpine/aports/-/issues/8626
pybind11_add_module(_prima NO_EXTRAS _prima.cpp)
target_include_directories(_prima PRIVATE ${CMAKE_SOURCE_DIR}/c/include)
target_link_libraries(_prima PRIVATE primac primaf)
target_compile_definitions(_prima PRIVATE VERSION_INFO=${PRIMA_VERSION})

# This section is cribbed from scipy/meson.build
if(WIN32 AND (CMAKE_C_COMPILER_ID MATCHES "GNU"))
  message(STATUS "MinGW detected, adding compiler and linker flags")
  target_link_options(_prima PUBLIC -lucrt -static)
  target_compile_options(_prima PUBLIC $<$<COMPILE_LANGUAGE:C>:-mlong-double-64 -D__USE_MINGW_ANSI_STDIO=1>)
  target_compile_options(_prima PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-D__USE_MINGW_ANSI_STDIO=1>)
endif()

install (TARGETS _prima DESTINATION prima/ COMPONENT Prima_Python_C_Extension)
# The following are not added to the component because scikit-build-core automatically
# detects files in python/<package_name> and adds them to the wheel. These commands
# are here in case one is building via cmake directly and not via scikit-build-core
file(GLOB SUPPORTING_PY_FILES "${CMAKE_CURRENT_SOURCE_DIR}/prima/*.py")
install (FILES ${SUPPORTING_PY_FILES} DESTINATION prima)

# We test the example files in CMake but not the test files. Those are tested by the Python
# build system (specifically cibuildwheel).
macro (prima_add_py_test name)
  add_test (NAME ${name}_example_python COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/examples/${name}_example.py)
  set_tests_properties (${name}_example_python PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_INSTALL_PREFIX}")
  # if (WIN32)
  #   file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}/bin" _BIN_PATH)
  #   set_property(TEST ${name}_example_python APPEND PROPERTY ENVIRONMENT "PATH=${_BIN_PATH}\\;$ENV{PATH}")
  # else ()
  #   set_property(TEST ${name}_example_python APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/c:$ENV{LD_LIBRARY_PATH};DYLD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/c:$ENV{DYLD_LIBRARY_PATH}")
  # endif ()
endmacro ()

prima_add_py_test (cobyla)
prima_add_py_test (bobyqa)
prima_add_py_test (newuoa)
prima_add_py_test (uobyqa)
prima_add_py_test (lincoa)
