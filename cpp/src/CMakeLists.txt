add_library (primacpp INTERFACE)
target_include_directories (primacpp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries (primacpp INTERFACE Eigen3::Eigen)
target_compile_options (primacpp INTERFACE $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)
target_compile_features (primacpp INTERFACE cxx_std_17)

include (GNUInstallDirs)
include (CMakePackageConfigHelpers)

install (
    TARGETS primacpp
    EXPORT primacpp-targets
)

install (
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/prima
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

install (
    EXPORT primacpp-targets
    FILE primacpp-targets.cmake
    NAMESPACE prima::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/primacpp
)

configure_package_config_file (
    ${CMAKE_CURRENT_SOURCE_DIR}/primacpp-config.cmake.in
    ${CMAKE_BINARY_DIR}/primacpp-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/primacpp
)

write_basic_package_version_file (
    ${CMAKE_BINARY_DIR}/primacpp-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install (
    FILES
        ${CMAKE_BINARY_DIR}/primacpp-config.cmake
        ${CMAKE_BINARY_DIR}/primacpp-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/primacpp
)
