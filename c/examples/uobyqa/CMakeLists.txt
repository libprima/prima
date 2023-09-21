
cmake_minimum_required(VERSION 3.13)
project(uobyqa_example C)
find_package(PRIMA CONFIG REQUIRED)
add_executable(uobyqa_example uobyqa_example.c)
target_link_libraries(uobyqa_example prima::primac)
install(TARGETS uobyqa_example DESTINATION bin)
