
cmake_minimum_required(VERSION 3.13)
project(cobyla_example C)
find_package(PRIMA CONFIG REQUIRED)
add_executable(cobyla_example cobyla_example.c)
target_link_libraries(cobyla_example prima::primac m)
install(TARGETS cobyla_example DESTINATION bin)
