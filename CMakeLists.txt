cmake_minimum_required(VERSION 3.22)

project(Darcy LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

if(NOT TARGET CommonDarcy)
  add_subdirectory(Common/Darcy)
endif()

ifem_add_application(
  NAME
    Darcy
  SOURCES
    DarcyArgs.C
    main_Darcy.C
  HEADERS
    DarcyArgs.h
  LIBRARIES
    CommonDarcy
)

enable_testing()

ifem_add_hdf5_test(
  TARGET
    Darcy
  TEST_FILES
    Tracer.hreg
    Wavefront_k10_p2_b20.hreg
)

ifem_add_hdf5_test(
  TARGET
    Darcy
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Lshape_p2_b5_adap.hreg
)

ifem_add_vtf_test(
  TARGET
    Darcy
  TEST_FILES
    Wavefront_k10_p2_b20.vreg
)

ifem_add_vtf_test(
  TARGET
    Darcy
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Lshape_p2_b5_adap.vreg
)

ifem_add_regression_test(
  TARGET
    Darcy
  CONDITIONS
    Zoltan_FOUND
  TEST_FILES
    MPI/Lshape_p2_b5.reg
  PARALLEL
    4
)

ifem_add_regression_test(
  TARGET
    Darcy
  CONDITIONS
    LRSpline_FOUND
    Zoltan_FOUND
  TEST_FILES
    MPI/Lshape_p2_b5_adap.reg
  PARALLEL
    4
)

ifem_add_regression_test(
  TARGET
    Darcy
  TEST_FILES
    Lshape_p2_b5.reg
    Wavefront_k10_p2_b20.reg
    Dirac.reg
    Dirac3D.reg
    Element.reg
    Element3D.reg
    SimpleSquare.reg
    SimpleSquare2p.reg
    Tracer.reg
    Tracer-be.reg
    Tracer-bdf2.reg
    Tracer-schedule.reg
)

ifem_add_regression_test(
  TARGET
    Darcy
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Lshape_p2_b5_adap.reg
)

foreach(adap pressure concentration total recovery recovery_press recovery_conc)
  ifem_add_regression_test(
    TARGET
      Darcy
    CONDITIONS
      LRSpline_FOUND
    TEST_FILES
      Tracer-adap-${adap}.reg
  )
endforeach()

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(TARGET Darcy)

# Installation
include(GNUInstallDirs)
install(TARGETS Darcy DESTINATION ${CMAKE_INSTALL_BINDIR})
