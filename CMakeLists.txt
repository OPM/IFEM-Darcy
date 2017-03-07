project(Darcy)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif()

# Generate regtest script with correct paths
configure_file(${IFEM_REGTEST_SCRIPT} regtest.sh)

set(DARCY_SRCS Darcy.C
               DarcyArgs.C
               DarcySolutions.C)

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR} ../Common)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_library(CommonDarcy STATIC ${DARCY_SRCS})

add_executable(Darcy main_Darcy.C)
target_link_libraries(Darcy CommonDarcy IFEMAppCommon ${IFEM_LIBRARIES})

# Installation
install(TARGETS Darcy DESTINATION bin COMPONENT bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")
add_doc_target(Darcy Darcy)

# Regression tests
enable_testing()
include(IFEMTesting)
if(MPI_FOUND)
  # Add parallel tests here
else()
  ifem_add_test(Lshape_p2_b5.reg Darcy)
  ifem_add_test(Wavefront_k10_p2_b20.reg Darcy)
endif()
list(APPEND TEST_APPS Darcy)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  Darcy
                  CommonDarcy IFEMAppCommon ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
