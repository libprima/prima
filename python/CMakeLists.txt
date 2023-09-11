
find_package (Python 3.6 COMPONENTS Interpreter)

execute_process(COMMAND ${Python_EXECUTABLE} -c "import scipy" RESULT_VARIABLE _IMPORT_SCIPY ERROR_QUIET)
if (NOT _IMPORT_SCIPY EQUAL 0)
  message(STATUS "scipy: NOT FOUND, not installing/testing Python bindings")
endif ()

if (Python_FOUND AND BUILD_SHARED_LIBS AND PRIMA_ENABLE_C AND (_IMPORT_SCIPY EQUAL 0))

  if (WIN32)
    set (PYTHON_SITE_PACKAGES Lib/site-packages)
  else ()
    set (PYTHON_SITE_PACKAGES ${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages)
  endif ()

  install (FILES prima.py DESTINATION ${PYTHON_SITE_PACKAGES})

  macro (prima_add_py_test name)
    add_test (NAME example_${name}_python COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/examples/example_${name}.py)
    set_tests_properties (example_${name}_python PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}")
    if (WIN32)
      file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}/bin" _BIN_PATH)
      set_property(TEST example_${name}_python APPEND PROPERTY ENVIRONMENT "PATH=${_BIN_PATH}\\;$ENV{PATH}")
    else ()
      set_property(TEST example_${name}_python APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/c:$ENV{LD_LIBRARY_PATH};DYLD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/c:$ENV{DYLD_LIBRARY_PATH}")
    endif ()
  endmacro ()

  prima_add_py_test (cobyla)
  prima_add_py_test (bobyqa)
  prima_add_py_test (newuoa)
  prima_add_py_test (uobyqa)
  prima_add_py_test (lincoa)

endif ()
