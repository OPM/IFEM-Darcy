project(Poisson)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

# Generate regtest script with correct paths
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)

set(DARCY_SRCS Darcy.C
               DarcySolutions.C)

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR} ../Common)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_library(CommonDarcy ${DARCY_SRCS})

add_executable(Darcy main_Darcy.C)
target_link_libraries(Darcy CommonDarcy ${IFEM_LIBRARIES})

# Installation
install(TARGETS Darcy DESTINATION bin COMPONENT bin)

enable_testing()

# Regression tests
if(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
else()
  file(GLOB DARCY_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                     ${PROJECT_SOURCE_DIR}/Test/*.reg)
  foreach(TESTFILE ${DARCY_TESTFILES})
    ifem_add_test(${TESTFILE} Darcy)
  endforeach()
endif()
list(APPEND TEST_APPS Darcy)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
