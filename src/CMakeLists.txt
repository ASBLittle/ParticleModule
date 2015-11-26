include_directories(${CMAKE_CURRENT_SOURCE_DIR})

project(vtkParticles)

set(FILE_SRCS vtkGmshReader.cxx vtkGranularTemperature.cxx vtkGenericProbe.cxx)
 
    
cmake_minimum_required(VERSION 2.8.8)

set(PACKAGE_VERSION 00.00.01)

FIND_PACKAGE(VTK REQUIRED)
FIND_PACKAGE(PythonLibs REQUIRED)
INCLUDE(${VTK_USE_FILE})
INCLUDE(${VTK_CMAKE_DIR}/vtkWrapPython.cmake)
find_program(PYTHON python)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(vtk_libs
    vtkCommon
    vtkFiltering
    vtkIO
    vtkGraphics)

set(VTK_MODULES_USED vtkCommon vtkFiltering vtkIO vtkGraphics)

link_directories(${VTK_DIR})
include_directories("${PYTHON_INCLUDE_DIRS}")

add_library(vtkParticles SHARED ${FILE_SRCS})
target_link_libraries(vtkParticles ${vtk_libs})

# Generate wrapper code
vtk_wrap_python3(vtkParticlesPython vtkParticlesPython_SRCS "${FILE_SRCS}")
    
add_library(vtkParticlesPythonD ${vtkParticlesPython_SRCS})
target_link_libraries(vtkParticlesPythonD  vtkParticles)   
    
#target_link_libraries(vtkParticlesPythonD         
#    ${VTK_LIBRARIES}
#    vtkPythonCore 
#    ${VTK_PYTHON_LIBRARIES} ${PYTHON_LIBRARY})
    
add_library(../particle_model/vtkParticlesPython MODULE ${vtkParticlesPython_SRCS} vtkParticlesPythonInit.cxx)
    
set(VTK_PYTHOND_LIBS)
foreach(TMP_LIB ${VTK_MODULES_USED})
     set(VTK_PYTHOND_LIBS ${VTK_PYTHOND_LIBS} ${TMP_LIB}PythonD)
endforeach()    
    
target_link_libraries(../particle_model/vtkParticlesPython vtkParticlesPythonD ${VTK_PYTHOND_LIBS})

set_target_properties(../particle_model/vtkParticlesPython PROPERTIES PREFIX "")

