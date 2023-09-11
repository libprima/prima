
cmake_minimum_required(VERSION 3.13)
project(lincoa_example C Fortran)
find_package(PRIMA CONFIG REQUIRED)
add_executable(lincoa_example lincoa_example.c)
target_link_libraries(lincoa_example prima::primac)
install(TARGETS lincoa_example DESTINATION bin)
