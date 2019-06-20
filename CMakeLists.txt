cmake_minimum_required(VERSION 2.8.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

PROJECT(particle_model)

set(CMAKE_MACOSX_RPATH 1)

set(USE_SYSTEM_VTK "ON" CACHE STRING "Use system VTK")
SET(Python_ADDITIONAL_VERSIONS "2.7 3.6 3.5 3.4")
set(PYTHON_VERSION 2.7 CACHE STRING "Python version to use to wrap VTK")
find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
FIND_PACKAGE(PythonLibs ${PYTHON_VERSION} REQUIRED)

if(USE_SYSTEM_VTK STREQUAL "ON")
  FIND_PACKAGE(VTK REQUIRED)  
else(USE_SYSTEM_VTK STREQUAL "ON")
  if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk/build/VTKTargets.cmake)
  else()
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vtk/build)
    set(VTK_CMAKE_ARGS "" CACHE STRING "Additional arguments to VTK build.")
    set(VTK_MAKE_JOBS "3" CACHE STRING "Number of jobs allowed for VTK build.")
    exec_program(${CMAKE_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/vtk/build 
      ARGS "-DVTK_PYTHON_VERSION=${PYTHON_VERSION} -DVTK_WRAP_PYTHON=ON -DVTK_USE_COCOA=OFF -DVTK_Group_MPI=ON -DBUILD_TESTING=OFF -DVTK_USE_X=ON ${VTK_CMAKE_ARGS} ${CMAKE_CURRENT_SOURCE_DIR}/vtk")
    exec_program(${CMAKE_MAKE_PROGRAM} ${CMAKE_CURRENT_BINARY_DIR}/vtk/build
      ARGS -j ${VTK_MAKE_JOBS})
  endif()
  set(VTK_DIR ${CMAKE_CURRENT_BINARY_DIR}/vtk/build)
endif(USE_SYSTEM_VTK STREQUAL "ON")

add_subdirectory(src)

enable_testing()
find_program(PYTEST_EXECUTABLE 
  NAMES pytest-${PYTHON_VERSION} py.test-${PYTHON_VERSION} py.test pytest)
add_test(NAME python_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${PYTEST_EXECUTABLE})
set_tests_properties(python_tests PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_CURRENT_BINARY_DIR}:particle_model FAIL_REGULAR_EXPRESSION "failed;error")

set(PACKAGE_VERSION 1.1.0)
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")	
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/particle_model/__init__.py")
set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

configure_file(${SETUP_PY_IN} ${SETUP_PY})

add_custom_command(OUTPUT ${OUTPUT}
                  COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build
                  COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                  DEPENDS ${DEPS})

add_custom_target(target ALL DEPENDS ${OUTPUT})

install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install)")
