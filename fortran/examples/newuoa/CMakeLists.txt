
cmake_minimum_required(VERSION 3.13)
project(newuoa_example Fortran)
find_package(PRIMA CONFIG REQUIRED)
add_executable(newuoa_example newuoa_example.f90)
target_link_libraries(newuoa_example prima::primaf)
install(TARGETS newuoa_example DESTINATION bin)
