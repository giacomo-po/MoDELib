cmake_minimum_required(VERSION 3.1.0)
#cmake_policy(SET CMP0057 NEW)
project(ddpy VERSION 1.0.0 LANGUAGES CXX)

option(USE_NOISE_GENERATOR "Use fftw to generate the noise" ON)

set(MoDELibRootDir "${CMAKE_CURRENT_SOURCE_DIR}/../..")
set(DDPY_INCLUDE_Dir "${CMAKE_CURRENT_SOURCE_DIR}")

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug) # this variable seems ineffective, so set -g manually
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")
set(CMAKE_CXX_FLAGS "-shared -std=gnu++2a -Wall -Wextra -march=native -fopenmp -Ofast -fPIC")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # creates compile_commands.json
set(CMAKE_VERBOSE_MAKEFILE TRUE)

execute_process(COMMAND "python3-config" "--extension-suffix" OUTPUT_VARIABLE PYTHON_SUFFIX)
#message("PYTHON_SUFFIX: ${PYTHON_SUFFIX}")

#set(DDPY "ddpy${PYTHON_SUFFIX}")
set(DDPY "ddpy")
#message("DDPY: ${DDPY}")
set(DDPY_OUTPUT_NAME "ddpy${PYTHON_SUFFIX}")
message("DDPY_OUTPUT_NAME IS: ${DDPY_OUTPUT_NAME}")
string(REPLACE ".so" ""  DDPY_OUTPUT_NAME ${DDPY_OUTPUT_NAME})
string(STRIP ${DDPY_OUTPUT_NAME} DDPY_OUTPUT_NAME)


# New FindPython mode  https://pybind11.readthedocs.io/en/stable/cmake/index.html
# either call find_package() or set PYBIND11_FINDPYTHON to ON
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
#pybind11_add_module(${DDPY} ddpy.cpp) # not working for cmake 3.16
Python_add_library(${DDPY} ddpy.cpp) # ???same source as above???
target_link_libraries(${DDPY} PUBLIC pybind11::headers Eigen3::Eigen pybind11::module ${DDOMP_Dir})
set_target_properties(${DDPY} PROPERTIES
                     CXX_VISIBILITY_PRESET default
                     VISIBILITY_INLINES_HIDDEN ON
                     LIBRARY_OUTPUT_NAME ${DDPY_OUTPUT_NAME}
                     PREFIX "")
   #INTERPROCEDURAL_OPTIMIZATION ON

#add_library(${DDPY} MODULE ddpy.cpp)
set_target_properties(${DDPY} PROPERTIES LIBRARY_OUTPUT_NAME ${DDPY_OUTPUT_NAME} PREFIX "")

# Classic mode  https://pybind11.readthedocs.io/en/stable/cmake/index.html
# set PythonLibsNew to influence python detection ???: set(PythonLibsNew TRUE) or put it in set_target_properties???
#find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
#find_package(pybind11 CONFIG REQUIRED)

message(STATUS "pybind11_FOUND: ${pybind11_FOUND}")

find_package(Eigen3 REQUIRED NAMES Eigen3)
if(USE_NOISE_GENERATOR)
	find_library(FFTW3_LIBRARIES fftw3)
	message(STATUS "FFTW3_INCLUDE_DIRS: ${FFTW3_INCLUDE_DIRS}")
	message(STATUS "FFTW3_LIBRARY_DIRS: ${FFTW3_LIBRARY_DIRS}")
	set(CMAKE_PREFIX_PATH ${FFTW3_LIBRARY_DIRS})
	include_directories(${FFTW3_INCLUDE_DIRS})
	include_directories(${BOOSTDir})
	add_definitions(-D _MODEL_GLIDE_PLANE_NOISE_GENERATOR_)
endif()

# append to DDPY_LINK_LIBRARIES_LIST: Eigen3::Eigen pybind11::module
list(APPEND DDPY_LINK_LIBRARIES_LIST Eigen3::Eigen pybind11::module)
if(USE_NOISE_GENERATOR)
   list(APPEND DDPY_LINK_LIBRARIES_LIST ${FFTW3_LIBRARIES})
endif()
message(STATUS "DDPY_LINK_LIBRARIES_LIST ${DDPY_LINK_LIBRARIES_LIST}")

target_link_libraries(${DDPY} PUBLIC ${DDPY_LINK_LIBRARIES_LIST})

target_include_directories(${DDPY} PUBLIC ${DDPY_INCLUDE_Dir})
target_include_directories(${DDPY} PUBLIC Eigen3::Eigen)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/Utilities/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/IO/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/Math/)
#target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/MPI/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/Geometry/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/Mesh/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/Lattices/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/FEM/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/GlidePlanes/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/LoopNetwork/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/Quadrature/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/ParticleInteraction/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DiscreteCrackMechanics/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/PolycrystallineMaterials/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DislocationDynamicsIO/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DislocationMicrostructure/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DislocationMobilities/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DislocationDynamicsBase/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DislocationDynamics/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DislocationDynamics/Polycrystals/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DislocationDynamics/Visualization/)
target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/DislocationDynamics/ElasticFields/)
#target_include_directories(${DDPY} PUBLIC ${MoDELibRootDir}/include/vtk/)

add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
add_subdirectory(${MoDELibRootDir}/src/GlidePlanes GlidePlanes)
add_subdirectory(${MoDELibRootDir}/src/DislocationMobilities DDDmobilities)
add_subdirectory(${MoDELibRootDir}/src/DislocationDynamicsIO DDDIO)
add_subdirectory(${MoDELibRootDir}/src/DislocationMicrostructure DDDmicrostructure)
add_subdirectory(${MoDELibRootDir}/src/PolycrystallineMaterials PolyMaterials)
add_subdirectory(${MoDELibRootDir}/src/DislocationDynamicsBase DDDBase)
add_subdirectory(${MoDELibRootDir}/src/DislocationDynamics DDD)
#add_subdirectory(${MoDELibRootDir}/src/vtk VTK)

