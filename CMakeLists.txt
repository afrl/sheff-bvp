cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0048 NEW)
project(FEAANALYSIS VERSION 1.0.0 LANGUAGES CXX)

# Options
option(FEAANALYSIS_BuildDev "Build the project in the dev folder. (intended for debugging this project by devs)" OFF)
option(FEAANALYSIS_BuildTests "Build the unit tests." OFF)

# For debugging, print generator, compiler, and options
message(STATUS "Gen: ${CMAKE_GENERATOR}, CXX: ${CMAKE_CXX_COMPILER_ID}, Mode: ${CMAKE_BUILD_TYPE}")
message(STATUS "Dev: ${FEAANALYSIS_BuildDev}, Tests: ${FEAANALYSIS_BuildTests}")

##################
# Checking input #
##################

include(${CMAKE_SOURCE_DIR}/cmake/CommonInputChecks.cmake)

########################
# Declare dependencies #
########################


find_package(Utility CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../utility/build" "${CMAKE_SOURCE_DIR}/../utility/install_dir")
find_package(SimpleMath CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../simple_math/build" "${CMAKE_SOURCE_DIR}/../simple_math/install_dir")
find_package(Constitutive CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../constitutive/build" "${CMAKE_SOURCE_DIR}/../constitutive/install_dir")
find_package(FlexMesh CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../flex_mesh/build" "${CMAKE_SOURCE_DIR}/../flex_mesh/install_dir")
find_package(FeaIo CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../fea_io/build" "${CMAKE_SOURCE_DIR}/../fea_io/install_dir")
find_package(FiniteElement CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../finite_element/build" "${CMAKE_SOURCE_DIR}/../finite_element/install_dir")
find_package(SparseLinearAlgebra CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../sparse_linear_algebra/build" "${CMAKE_SOURCE_DIR}/../sparse_linear_algebra/install_dir")
find_package(Output CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../output/build" "${CMAKE_SOURCE_DIR}/../output/install_dir")
find_package(DofManager CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../dof_manager/build" "${CMAKE_SOURCE_DIR}/../dof_manager/install_dir")
find_package(Assembler CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../assembler/build" "${CMAKE_SOURCE_DIR}/../assembler/install_dir")
find_package(FeaCrack CONFIG REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../fea_crack/build" "${CMAKE_SOURCE_DIR}/../fea_crack/install_dir")

################################
# Setting up FEAANALYSIS target #
################################

# Set the variables for FEAANALYSIS
set(FEAANALYSIS_INCLUDE ${CMAKE_SOURCE_DIR}/include)
set(FEAANALYSIS_INCLUDE_SUB ${FEAANALYSIS_INCLUDE}/fea_analysis)
set(FEAANALYSIS_SRC ${CMAKE_SOURCE_DIR}/source)
file (GLOB_RECURSE FEAANALYSIS_HEADER_LIST "${FEAANALYSIS_INCLUDE}/*.h")
file (GLOB_RECURSE FEAANALYSIS_SOURCE_LIST "${FEAANALYSIS_SRC}/*.h" "${FEAANALYSIS_SRC}/*.cpp")
    
# Create target
add_library(feaanalysis STATIC)
# Add the header includes and headers to the library
target_include_directories(feaanalysis PUBLIC
    $<BUILD_INTERFACE:${FEAANALYSIS_INCLUDE}>
    $<INSTALL_INTERFACE:include>)
# Add headers in the include folder as PUBLIC so downstream
# targets have access to them
target_sources(feaanalysis PUBLIC ${FEAANALYSIS_HEADER_LIST})
# Add headers and sources in the source folder as PRIVATE as they
# will not be needed downstream
target_sources(feaanalysis PRIVATE ${FEAANALYSIS_SOURCE_LIST})
# Add an alias so that a library can be used inside the build tree, eg when testing
add_library(FeaAnalysis::FeaAnalysis ALIAS feaanalysis)
# Add C++ 14 flag
target_compile_features(feaanalysis PUBLIC cxx_std_14)
# Link dependencies
target_link_libraries(feaanalysis PUBLIC Utility::Utility)
target_link_libraries(feaanalysis PUBLIC SimpleMath::SimpleMath)
target_link_libraries(feaanalysis PUBLIC Constitutive::Constitutive)
target_link_libraries(feaanalysis PUBLIC FlexMesh::FlexMesh)
target_link_libraries(feaanalysis PUBLIC FeaIo::FeaIo)
target_link_libraries(feaanalysis PUBLIC FiniteElement::FiniteElement)
target_link_libraries(feaanalysis PUBLIC SparseLinearAlgebra::SparseLinearAlgebra)
target_link_libraries(feaanalysis PUBLIC Output::Output)
target_link_libraries(feaanalysis PUBLIC DofManager::DofManager)
target_link_libraries(feaanalysis PUBLIC Assembler::Assembler)
target_link_libraries(feaanalysis PUBLIC FeaCrack::FeaCrack)


# Assign source groups
include(${CMAKE_SOURCE_DIR}/cmake/Utilities.cmake)
assign_source_group("fea_analysis" ${FEAANALYSIS_HEADER_LIST})
assign_source_group("fea_analysis" ${FEAANALYSIS_SOURCE_LIST})
assign_source_group_Utility()
assign_source_group_SimpleMath()
assign_source_group_Constitutive()
assign_source_group_FlexMesh()
assign_source_group_FeaIo()
assign_source_group_FiniteElement()
assign_source_group_SparseLinearAlgebra()
assign_source_group_Output()
assign_source_group_DofManager()
assign_source_group_Assembler()
assign_source_group_FeaCrack()


#################################
# Add dev and unit test targets #
#################################

# Build dev
if(FEAANALYSIS_BuildDev)
    add_subdirectory(dev)
endif()
    
# Build tests
if(FEAANALYSIS_BuildTests)
    include(CTest)
    enable_testing()
    add_subdirectory(unit_tests)
endif()

################
# Installation #
################

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/FeaAnalysis)

# Install the feaanalysis target (library and headers)
install(
    TARGETS feaanalysis
    EXPORT FeaAnalysisTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)
# This is required so that the exported target has the name FeaAnalysis and not feaanalysis
set_target_properties(feaanalysis PROPERTIES EXPORT_NAME FeaAnalysis)

# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/FeaAnalysisConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Create the package config file for find_package(FeaAnalysis)
configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FeaAnalysisConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/FeaAnalysisConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Install the config and config version files
# Note: if you write custom find modules, add them here
# Refer to https://github.com/pabloariasal/modern-cmake-sample/blob/cfda7b51ff915409fd61d4f2e75b63032bedb182/libjsonutils/CMakeLists.txt
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/FeaAnalysisConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/FeaAnalysisConfigVersion.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/FeaAnalysisTargets.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

###########
# Exports #
###########

# Export just the feaanalysis target (not the dev, tests, etc.)
export(
    TARGETS feaanalysis
    FILE ${CMAKE_CURRENT_BINARY_DIR}/FeaAnalysisTargets.cmake
    NAMESPACE FeaAnalysis::
)

# Register package in user's package registry
export(PACKAGE FeaAnalysis)
