cmake_minimum_required(VERSION 3.22)

project(CommonDarcy LANGUAGES C CXX)

ifem_add_library(
  NAME
    CommonDarcy
  SOURCES
    Darcy.C
    DarcyAdvection.C
    DarcyMaterial.C
    DarcySolutions.C
    DarcyTransport.C
    SIMDarcy.C
    SIMDarcyAdvection.C
    SIMDarcySchedule.C
  HEADERS
    Darcy.h
    DarcyAdvection.h
    DarcyEnums.h
    DarcyMaterial.h
    DarcySolutions.h
    DarcyTransport.h
    SIMDarcy.h
    SIMDarcyAdap.h
    SIMDarcyAdvection.h
    SIMDarcySchedule.h
  LIBRARIES
    IFEMAppCommon
    IFEM
)

# Unit tests
ifem_add_test_app(
  SOURCES
    ${PROJECT_SOURCE_DIR}/Test/TestSIMDarcy.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/../../Test
  NAME
    Darcy
  LIBRARIES
    CommonDarcy
)

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
