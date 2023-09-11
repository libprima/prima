
cmake_minimum_required(VERSION 3.13)
project(bobyqa_example Fortran)
find_package(PRIMA CONFIG REQUIRED)
add_executable(bobyqa_example bobyqa_example.f90)
target_link_libraries(bobyqa_example prima::primaf)
install(TARGETS bobyqa_example DESTINATION bin)
