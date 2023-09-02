cmake_minimum_required(VERSION 3.5)

project(Darcy)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif()

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR} ../Common Common/Darcy)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

enable_testing()
include(IFEMTesting)

add_subdirectory(Common/Darcy)

add_executable(Darcy main_Darcy.C DarcyArgs.C)
target_link_libraries(Darcy CommonDarcy IFEMAppCommon ${IFEM_LIBRARIES})
list(APPEND CHECK_SOURCES ${DARCY_SRCS} main_Darcy.C DarcyArgs.C)

# Installation
install(TARGETS Darcy DESTINATION bin COMPONENT bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")
add_doc_target(Darcy Darcy)

# Regression tests
if(MPI_FOUND)
  ifem_add_test(MPI/Lshape_p2_b5.reg Darcy 4)
  if(LRSPLINE_FOUND OR LRSpline_FOUND)
    ifem_add_test(MPI/Lshape_p2_b5_adap.reg Darcy 4)
  endif()
endif()
if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  ifem_add_test(Lshape_p2_b5.reg Darcy)
  ifem_add_test(Wavefront_k10_p2_b20.reg Darcy)
  ifem_add_test(Dirac.reg Darcy)
  ifem_add_test(Dirac3D.reg Darcy)
  ifem_add_test(Element.reg Darcy)
  ifem_add_test(Element3D.reg Darcy)
  ifem_add_test(Tracer.reg Darcy)
  ifem_add_test(Tracer-be.reg Darcy)
  ifem_add_test(Tracer-bdf2.reg Darcy)
  ifem_add_test(Tracer-schedule.reg Darcy)
  if(LRSPLINE_FOUND OR LRSpline_FOUND)
    foreach(adap pressure concentration total recovery recovery_press recovery_conc)
      ifem_add_test(Tracer-adap-${adap}.reg Darcy)
    endforeach()
  endif()
  ifem_add_test(SimpleSquare.reg Darcy)
  ifem_add_test(SimpleSquare2p.reg Darcy)
  ifem_add_vtf_test(Wavefront_k10_p2_b20.vreg Darcy)
  ifem_add_hdf5_test(Wavefront_k10_p2_b20.hreg Darcy)
  if(LRSPLINE_FOUND OR LRSpline_FOUND)
    ifem_add_test(Lshape_p2_b5_adap.reg Darcy)
    ifem_add_hdf5_test(Lshape_p2_b5_adap.hreg Darcy)
    ifem_add_vtf_test(Lshape_p2_b5_adap.vreg Darcy)
  endif()
  ifem_add_hdf5_test(Tracer.hreg Darcy)
endif()
list(APPEND TEST_APPS Darcy)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  add_check_target()
endif()
