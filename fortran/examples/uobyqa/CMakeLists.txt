
cmake_minimum_required(VERSION 3.13)
project(uobyqa_example Fortran)
find_package(PRIMA CONFIG REQUIRED)
add_executable(uobyqa_example uobyqa_example.f90)
target_link_libraries(uobyqa_example prima::primaf)
install(TARGETS uobyqa_example DESTINATION bin)
