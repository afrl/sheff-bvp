cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0048 NEW)
project(BVP VERSION 1.0.0 LANGUAGES CXX)

# Set package toolchain location for Linux. To enable efficient switching between
# toolchains for cross-compiling, we modify the location where CMake looks for 
# packages using the find_package command by symbolically linking the default location,
# $HOME/.cmake/packages, to a similar location with the folder name given by the 
# compiler family. 
if (UNIX)
    execute_process(COMMAND "mkdir" "-p" "$ENV{HOME}/.cmake/${CMAKE_CXX_COMPILER_ID}"
      RESULT_VARIABLE toolchain_res OUTPUT_VARIABLE toolchain_out)
    if (toolchain_res)
        message(ERROR "Couldn't create toolchain dir: $ENV{HOME}/.cmake/${CMAKE_CXX_COMPILER_ID} ${toolchain_out}")
    endif ()
    execute_process(COMMAND "ln" "-sfn" "$ENV{HOME}/.cmake/${CMAKE_CXX_COMPILER_ID}" "$ENV{HOME}/.cmake/packages"
      RESULT_VARIABLE toolchain_res OUTPUT_VARIABLE toolchain_out)
    if (toolchain_res)
        message(ERROR "Couldn't create toolchain link: $ENV{HOME}/.cmake/packages ${toolchain_out}")
    endif ()
endif()


# Options
option(BVP_BuildTests "Build the unit tests." OFF)

# For debugging, print generator, compiler, and options
message(STATUS "Gen: ${CMAKE_GENERATOR}, CXX: ${CMAKE_CXX_COMPILER_ID}, Mode: ${CMAKE_BUILD_TYPE}")
message(STATUS "Tests: ${BVP_BuildTests}")

# Add cmake dir to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

##################
# Checking input #
##################

include(${CMAKE_SOURCE_DIR}/cmake/CommonInputChecks.cmake)

########################
# Declare dependencies #
########################


find_package(Utility CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../utility/install_dir" "${CMAKE_SOURCE_DIR}/../../utility/install_dir" "${CMAKE_SOURCE_DIR}/modules/utility/install_dir")
find_package(SimpleMath CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../simple_math/install_dir" "${CMAKE_SOURCE_DIR}/../../simple_math/install_dir" "${CMAKE_SOURCE_DIR}/modules/simple_math/install_dir")
find_package(Constitutive CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../constitutive/install_dir" "${CMAKE_SOURCE_DIR}/../../constitutive/install_dir" "${CMAKE_SOURCE_DIR}/modules/constitutive/install_dir")
  find_package(Mesh CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../mesh/install_dir" "${CMAKE_SOURCE_DIR}/../../mesh/install_dir" "${CMAKE_SOURCE_DIR}/modules/mesh/install_dir")
find_package(FeaIo CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../fea_io/install_dir" "${CMAKE_SOURCE_DIR}/../../fea_io/install_dir" "${CMAKE_SOURCE_DIR}/modules/fea_io/install_dir")
find_package(FiniteElement CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../finite_element/install_dir" "${CMAKE_SOURCE_DIR}/../../finite_element/install_dir" "${CMAKE_SOURCE_DIR}/modules/finite_element/install_dir")
find_package(SparseLinearAlgebra CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../sparse_linear_algebra/install_dir" "${CMAKE_SOURCE_DIR}/../../sparse_linear_algebra/install_dir" "${CMAKE_SOURCE_DIR}/modules/sparse_linear_algebra/install_dir")
find_package(DofManager CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../dof_manager/install_dir" "${CMAKE_SOURCE_DIR}/../../dof_manager/install_dir" "${CMAKE_SOURCE_DIR}/modules/dof_manager/install_dir")
find_package(Assembler CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../assembler/install_dir" "${CMAKE_SOURCE_DIR}/../../assembler/install_dir" "${CMAKE_SOURCE_DIR}/modules/assembler/install_dir")
find_package(Xfem CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../xfem/install_dir" "${CMAKE_SOURCE_DIR}/../../xfem/install_dir" "${CMAKE_SOURCE_DIR}/modules/xfem/install_dir")


################################
# Setting up BVP target #
################################

# Set the variables for BVP
set(BVP_INCLUDE ${CMAKE_SOURCE_DIR}/include)
set(BVP_SRC ${CMAKE_SOURCE_DIR}/source)
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


# Assign source groups
include(${CMAKE_SOURCE_DIR}/cmake/Utilities.cmake)
assign_source_group("bvp" ${BVP_HEADER_LIST})
assign_source_group("bvp" ${BVP_SOURCE_LIST})
assign_source_group_Utility()
assign_source_group_SimpleMath()
assign_source_group_Constitutive()
assign_source_group_Mesh()
assign_source_group_FeaIo()
assign_source_group_FiniteElement()
assign_source_group_SparseLinearAlgebra()
assign_source_group_DofManager()
assign_source_group_Assembler()
assign_source_group_Xfem()
assign_source_group_SecondPartyMath()
assign_source_group_SecondPartyMesh()
assign_source_group_SecondPartyMetis()

#################################
# Add unit test targets #
#################################

# Build tests
if(BVP_BuildTests)
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

# Register package in user's package registry
export(PACKAGE Bvp)
