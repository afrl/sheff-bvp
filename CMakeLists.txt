cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0048 NEW)
project(BVP VERSION 1.0.0 LANGUAGES CXX)

# Options
option(SHEFF_BuildTests "Build the unit tests." OFF)

# For debugging, print generator, compiler, and options
message(STATUS "Gen: ${CMAKE_GENERATOR}, CXX: ${CMAKE_CXX_COMPILER_ID}, Mode: ${CMAKE_BUILD_TYPE}")
message(STATUS "Tests: ${SHEFF_BuildTests}")

# Add cmake dir to module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

##################
# Checking input #
##################

include(${PROJECT_SOURCE_DIR}/cmake/CommonInputChecks.cmake)

########################
# Declare dependencies #
########################

if(NOT TARGET Assembler::Assembler)
  message(STATUS "Assembler did not already exist, finding package.")
  find_package(Assembler CONFIG REQUIRED NO_CMAKE_SYSTEM_PACKAGE_REGISTRY HINTS
               "${PROJECT_SOURCE_DIR}/../../first_party/assembler/build")
endif()

if(NOT TARGET Xfem::Xfem)
  message(STATUS "Xfem did not already exist, finding package.")
  find_package(Xfem CONFIG REQUIRED NO_CMAKE_SYSTEM_PACKAGE_REGISTRY HINTS
               "${PROJECT_SOURCE_DIR}/../../first_party/xfem/build")
endif()

################################
# Setting up BVP target #
################################

# Set the variables for BVP
set(BVP_INCLUDE ${PROJECT_SOURCE_DIR}/include)
set(BVP_SRC ${PROJECT_SOURCE_DIR}/source)
file (GLOB_RECURSE BVP_HEADER_LIST "${BVP_INCLUDE}/*.h")
file (GLOB_RECURSE BVP_SOURCE_LIST "${BVP_SRC}/*.h" "${BVP_SRC}/*.cpp")
    
# Create target
add_library(bvp STATIC)
# Add the header includes and headers to the library
target_include_directories(bvp PUBLIC
    $<BUILD_INTERFACE:${BVP_INCLUDE}>
    $<INSTALL_INTERFACE:include>)
# Add headers in the include folder as PUBLIC so downstream
# targets have access to them
target_sources(bvp PUBLIC ${BVP_HEADER_LIST})
# Add headers and sources in the source folder as PRIVATE as they
# will not be needed downstream
target_sources(bvp PRIVATE ${BVP_SOURCE_LIST})
# Add an alias so that a library can be used inside the build tree, eg when testing
add_library(Bvp::Bvp ALIAS bvp)
# Add C++ 14 flag
target_compile_features(bvp PUBLIC cxx_std_14)
set_target_properties(bvp PROPERTIES FOLDER "sheff")
# Link dependencies
target_link_libraries(bvp PUBLIC Utility::Utility)
target_link_libraries(bvp PUBLIC SimpleMath::SimpleMath)
target_link_libraries(bvp PUBLIC Constitutive::Constitutive)
target_link_libraries(bvp PUBLIC Mesh::Mesh)
target_link_libraries(bvp PUBLIC FeaIo::FeaIo)
target_link_libraries(bvp PUBLIC FiniteElement::FiniteElement)
target_link_libraries(bvp PUBLIC SparseLinearAlgebra::SparseLinearAlgebra)
target_link_libraries(bvp PUBLIC DofManager::DofManager)
target_link_libraries(bvp PUBLIC Assembler::Assembler)
target_link_libraries(bvp PUBLIC Xfem::Xfem)

###########################
# Assigning source groups #
###########################

include("$ENV{SHEFF_TOOLS}/cmake/IncludeSHEFF.cmake")
register_sources_sheff("bvp"
                       "${BVP_HEADER_LIST}"
                       "${BVP_SOURCE_LIST}")
assign_source_groups_sheff()

#################################
# Add unit test targets #
#################################

# Build tests
if(SHEFF_BuildTests)
    include(CTest)
    enable_testing()
    add_subdirectory(unit_tests)
endif()

################
# Installation #
################

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Bvp)

# Install the bvp target (library and headers)
install(
    TARGETS bvp
    EXPORT BvpTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)
# This is required so that the exported target has the name Bvp and not bvp
set_target_properties(bvp PROPERTIES EXPORT_NAME Bvp)

# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/BvpConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Create the package config file for find_package(Bvp)
configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/BvpConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/BvpConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Install the config and config version files
# Note: if you write custom find modules, add them here
# Refer to https://github.com/pabloariasal/modern-cmake-sample/blob/cfda7b51ff915409fd61d4f2e75b63032bedb182/libjsonutils/CMakeLists.txt
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/BvpConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/BvpConfigVersion.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/BvpTargets.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

###########
# Exports #
###########

# Export just the bvp target (not tests, etc.)
export(
    TARGETS bvp
    FILE ${CMAKE_CURRENT_BINARY_DIR}/BvpTargets.cmake
    NAMESPACE Bvp::
)
