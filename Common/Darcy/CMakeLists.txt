cmake_minimum_required(VERSION 3.5)

project(CommonDarcy)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

include_directories(../../../Common)

set(DARCY_SOURCES Darcy.C
                  Darcy.h
                  DarcyAdvection.C
                  DarcyAdvection.h
                  DarcyEnums.h
                  DarcyMaterial.C
                  DarcyMaterial.h
                  DarcySolutions.C
                  DarcySolutions.h
                  DarcyTransport.C
                  DarcyTransport.h
                  SIMDarcy.C
                  SIMDarcy.h
                  SIMDarcyAdap.h
                  SIMDarcyAdvection.C
                  SIMDarcyAdvection.h
                  SIMDarcySchedule.C
                  SIMDarcySchedule.h)

add_library(CommonDarcy STATIC ${DARCY_SOURCES})

# Unit tests
include(IFEMTesting)

IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/../../Test
                  Darcy 0
                  CommonDarcy IFEMAppCommon ${IFEM_LIBRARIES})


foreach(src ${DARCY_SOURCES})
  list(APPEND CHECK_SOURCES ${PROJECT_SOURCE_DIR}/${src})
endforeach()

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
