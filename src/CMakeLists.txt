include_directories(${CMAKE_CURRENT_SOURCE_DIR})

project(vtkParticles)

set(FILE_SRCS vtkGmshReader.cxx vtkGranularTemperature.cxx vtkGenericProbe.cxx vtkShowCVs.cxx vtkMergePointFilter.cxx)
    
cmake_minimum_required(VERSION 2.8.8)

set(PACKAGE_VERSION 00.00.01)

set(VTK_LIBRARY_SUFFIX CACHE STRING "VTK library suffix")
FIND_PACKAGE(VTK REQUIRED)
SET(Python_ADDITIONAL_VERSIONS "2.7 3.6 3.5 3.4")
set(PYTHON_VERSION 2.7 CACHE STRING "Python version to use")
find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
FIND_PACKAGE(PythonLibs ${PYTHON_VERSION} REQUIRED)
INCLUDE(${VTK_USE_FILE})
INCLUDE(${VTK_CMAKE_DIR}/vtkWrapPython.cmake)
find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import imp\ntry:\n print(imp.find_module('numpy')[1]+'/core/include/')\nexcept:\n pass"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR)
if (NUMPY_INCLUDE_DIR STREQUAL "")
  set(NUMPY_INCLUDE_DIR NUMPY_INCLUDE_DIR-NOTFOUND)
  find_path(NUMPY_INCLUDE_DIR numpy/arrayobject.h)
endif (NUMPY_INCLUDE_DIR STREQUAL "")
include_directories(${NUMPY_INCLUDE_DIR})

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

if (VTK_MAJOR_VERSION LESS 6)
  set(VTK_MODULES_USED vtkFiltering vtkIO vtkGraphics)
  set(PYTHOND_SUFFIX PythonD)
else(VTK_MAJOR_VERSION LESS 6)
  set(VTK_MODULES_USED vtkFiltersCore vtkCommonExecutionModel vtkCommonCore)
  set(PYTHOND_SUFFIX Python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}D${VTK_LIBRARY_SUFFIX}-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION})
endif()

link_directories(${VTK_DIR})
link_directories(${VTK_DIR}/lib)
link_directories(/opt/local/lib)
include_directories("${PYTHON_INCLUDE_DIRS}")

execute_process(
  COMMAND git log -1 --format=%H
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION_STRING
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

add_definitions(-DGIT_VERSION_STRING="${GIT_VERSION_STRING}")

add_library(vtkParticles SHARED ${FILE_SRCS})
target_link_libraries(vtkParticles ${VTK_LIBRARIES})
add_library(vtk_extras SHARED Picker.cxx BoundingSurface.cxx)
target_link_libraries(vtk_extras ${VTK_LIBRARIES})

# Generate wrapper code
vtk_wrap_python3(vtkParticlesPython vtkParticlesPython_SRCS "${FILE_SRCS}")
    
add_library(vtkParticlesPythonD ${vtkParticlesPython_SRCS})
target_link_libraries(vtkParticlesPythonD  vtkParticles)   
    
target_link_libraries(vtkParticlesPythonD         
    ${VTK_LIBRARIES}
    ${VTK_PYTHON_LIBRARIES} ${PYTHON_LIBRARY})
    
add_library(vtkParticlesPython MODULE ${vtkParticlesPython_SRCS} vtkParticlesPythonInit.cxx)

add_library(vtk_extras_mod MODULE PickerObject.cxx vtk_extras.cxx)
        
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(DEFINED VTK_CUSTOM_LIBRARY_SUFFIX)
  set(VTK_LIB_SUFFIX ${VTK_CUSTOM_LIBRARY_SUFFIX})
else()
  set(VTK_LIB_SUFFIX -${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION})
endif()

set(VTK_PYTHOND_LIBS)
foreach(TMP_LIB ${VTK_MODULES_USED})
     set(VTK_PYTHOND_LIBS ${VTK_PYTHOND_LIBS} ${TMP_LIB}${PYTHOND_SUFFIX})
endforeach()    

target_link_libraries(vtkParticlesPython vtkParticlesPythonD ${VTK_PYTHON_LIBRARIES} ${VTK_PYTHOND_LIBS})
target_link_libraries(vtk_extras_mod vtkParticlesPythonD ${VTK_PYTHON_LIBRARIES} ${VTK_PYTHOND_LIBS} vtk_extras)

set_target_properties(vtkParticlesPython PROPERTIES PREFIX "")
set_target_properties(vtk_extras_mod PROPERTIES PREFIX "")

ADD_CUSTOM_COMMAND(TARGET vtk_extras_mod POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy  ${CMAKE_CURRENT_BINARY_DIR}/vtk_extras_mod.so ${CMAKE_CURRENT_BINARY_DIR}/../particle_model/vtk_extras.so)
ADD_CUSTOM_COMMAND(TARGET vtkParticlesPython POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy  ${CMAKE_CURRENT_BINARY_DIR}/vtkParticlesPython.so ${CMAKE_CURRENT_BINARY_DIR}/../particle_model/vtkParticlesPython.so)

