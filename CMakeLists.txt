cmake_minimum_required (VERSION 3.13)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (prima C Fortran)

option (BUILD_SHARED_LIBS "shared/static" ON)
if (WIN32)
  # symbols are not exported on windows
  set (BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
endif ()

include (GNUInstallDirs)

# the compiler enables an executable stack because of nested functions and this is fine
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  include (CheckLinkerFlag)
  check_linker_flag (Fortran "-Wl,--no-warn-execstack" HAVE_WARN_EXECSTACK)
endif ()

# allocate arrays on the heap instead of the stack, slower but can avoid memory errors on large problems
option (PRIMA_HEAP_ARRAYS "allocate arrays on heap" ON)
if (PRIMA_HEAP_ARRAYS)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-stack-arrays")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if (WIN32)
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /heap-arrays")
    else ()
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -heap-arrays")
    endif ()
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mnostack_arrays")
  endif ()
endif ()

option(PRIMA_ENABLE_EXAMPLES "build examples by default" OFF)
add_custom_target (examples)
enable_testing ()

option(PRIMA_ENABLE_TESTING "build tests" OFF)
add_custom_target (tests)
add_dependencies(tests examples)
add_subdirectory(fortran)

option (PRIMA_ENABLE_C "C binding" ON)

if (PRIMA_ENABLE_C)
  add_subdirectory(c)
  set(primac_target "primac")
endif ()

file (READ VERSION.txt Prima_VERSION)


install(
    TARGETS primaf ${primac_target}
    EXPORT prima-targets
    INCLUDES DESTINATION include
)

install(
    EXPORT prima-targets
    FILE prima-targets.cmake
    NAMESPACE prima::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prima
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/prima-config.cmake.in
    ${CMAKE_BINARY_DIR}/prima-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prima
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/prima-config-version.cmake
    VERSION ${Prima_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/prima-config.cmake
        ${CMAKE_BINARY_DIR}/prima-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prima
)
