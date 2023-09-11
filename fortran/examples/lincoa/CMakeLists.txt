
cmake_minimum_required(VERSION 3.13)
project(lincoa_example Fortran)
find_package(PRIMA CONFIG REQUIRED)
add_executable(lincoa_example lincoa_example.f90)
target_link_libraries(lincoa_example prima::primaf)
install(TARGETS lincoa_example DESTINATION bin)
