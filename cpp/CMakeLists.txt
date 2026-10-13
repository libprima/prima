cmake_minimum_required (VERSION 3.18)

project (primacpp VERSION 0.1.0 LANGUAGES CXX)

option (BUILD_SHARED_LIBS "shared/static" ON)

include (GNUInstallDirs)

find_package (Eigen3 REQUIRED)
message (STATUS "Found Eigen3: ${Eigen3_DIR} (found version ${Eigen3_VERSION})")

add_subdirectory (src)
enable_testing()
add_subdirectory (tests)
add_subdirectory (examples)
