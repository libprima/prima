
cmake_minimum_required(VERSION 3.13)
project(newuoa_example C Fortran)
find_package(PRIMA CONFIG REQUIRED)
add_executable(newuoa_example newuoa_example.c)
target_link_libraries(newuoa_example prima::primac)
install(TARGETS newuoa_example DESTINATION bin)
