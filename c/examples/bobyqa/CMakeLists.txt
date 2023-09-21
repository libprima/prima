
cmake_minimum_required(VERSION 3.13)
project(bobyqa_example C)
find_package(PRIMA CONFIG REQUIRED)
add_executable(bobyqa_example bobyqa_example.c)
target_link_libraries(bobyqa_example prima::primac)
install(TARGETS bobyqa_example DESTINATION bin)
