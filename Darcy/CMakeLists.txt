cmake_minimum_required(VERSION 3.22)

project(DarcyLib LANGUAGES C CXX)

ifem_add_library(
  NAME
    Drcy
  ALIAS
    IFEM::Darcy
  SOURCES
    Darcy.C
    DarcyAdvection.C
    DarcyMaterial.C
    DarcySolutions.C
    DarcyTransport.C
    DarcyTransportCorr.C
    SIMDarcy.C
    SIMDarcyAdvection.C
    SIMDarcySchedule.C
    SIMDarcyTransportCorr.C
  HEADERS
    Darcy.h
    DarcyAdvection.h
    DarcyEnums.h
    DarcyMaterial.h
    DarcySolutions.h
    DarcyTransport.h
    DarcyTransportCorr.h
    SIMDarcy.h
    SIMDarcyAdap.h
    SIMDarcyAdvection.h
    SIMDarcySchedule.h
    SIMDarcyTransportCorr.h
  LIBRARIES
    IFEM
    IFEM::AppUtils
)
set_target_properties(Drcy
  PROPERTIES
  ARCHIVE_OUTPUT_NAME
    Darcy
)

# Unit tests
ifem_add_test_app(
  NAME
    Darcy
  SOURCES
    ${PROJECT_SOURCE_DIR}/Test/TestSIMDarcy.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/../Test
  LIBRARIES
    IFEM::Darcy
)

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
