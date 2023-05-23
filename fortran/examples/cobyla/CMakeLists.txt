
cmake_minimum_required(VERSION 3.13)
project(cobyla_example Fortran)
find_package(PRIMA CONFIG REQUIRED)
add_executable(cobyla_example cobyla_example.f90)
target_link_libraries(cobyla_example prima::primaf)
install(TARGETS cobyla_example DESTINATION bin)
