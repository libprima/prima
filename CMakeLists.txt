cmake_minimum_required (VERSION 3.13)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (prima Fortran)

option (BUILD_SHARED_LIBS "shared/static" ON)

include (GNUInstallDirs)

# the compiler enables an executable stack because of nested functions and this is fine
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  include (CheckLinkerFlag)
  check_linker_flag (Fortran "-Wl,--no-warn-execstack" HAVE_WARN_EXECSTACK)
endif ()

# Set additional Fortran compiler flags
# 1. We require the compilers to allocate arrays on the heap instead of the stack, which is
# slower (does not matter for DFO applications) but can avoid memory errors on large problems.
# 2. We require the compilers to compile the solvers so that they can be called recursively.
# See https://fortran-lang.discourse.group/t/frecursive-assume-recursion-and-recursion-thread-safety
option (PRIMA_HEAP_ARRAYS "allocate arrays on heap" ON)
if (PRIMA_HEAP_ARRAYS)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-stack-arrays -frecursive")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if (WIN32)
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /heap-arrays /assume:recursion")
    else ()
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -heap-arrays -assume recursion")
    endif ()
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mnostack_arrays -Mrecursive")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive")  # What about stack/heap?
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
  enable_language(C)
  add_subdirectory(c)
  set(primac_target "primac")
endif ()

file (STRINGS VERSION.txt PRIMA_VERSION)

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
    VERSION ${PRIMA_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/prima-config.cmake
        ${CMAKE_BINARY_DIR}/prima-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prima
)
