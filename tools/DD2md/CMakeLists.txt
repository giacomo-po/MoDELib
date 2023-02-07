cmake_minimum_required(VERSION 3.1.0)
#cmake_policy(SET CMP0057 NEW)
project(dd2md VERSION 1.0.0 LANGUAGES CXX)

option(USE_NOISE_GENERATOR "Use fftw to generate the noise" ON)

set(MoDELibRootDir "${CMAKE_CURRENT_SOURCE_DIR}/../..")
set(DD2MD_INCLUDE_Dir "${CMAKE_CURRENT_SOURCE_DIR}")

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug) # this variable seems ineffective, so set -g manually
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")
set(CMAKE_CXX_FLAGS "-shared -std=gnu++2a -Wall -Wextra -march=native -fopenmp -Ofast -fPIC -g")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # creates compile_commands.json
set(CMAKE_VERBOSE_MAKEFILE TRUE)

execute_process(COMMAND "python3-config" "--extension-suffix" OUTPUT_VARIABLE PYTHON_SUFFIX)
#message("PYTHON_SUFFIX: ${PYTHON_SUFFIX}")

#set(DD2MD "dd2md${PYTHON_SUFFIX}")
set(DD2MD "dd2md")
#message("DD2MD: ${DD2MD}")
set(DD2MD_OUTPUT_NAME "dd2md${PYTHON_SUFFIX}")
message("DD2MD_OUTPUT_NAME IS: ${DD2MD_OUTPUT_NAME}")
string(REPLACE ".so" ""  DD2MD_OUTPUT_NAME ${DD2MD_OUTPUT_NAME})
string(STRIP ${DD2MD_OUTPUT_NAME} DD2MD_OUTPUT_NAME)


# New FindPython mode  https://pybind11.readthedocs.io/en/stable/cmake/index.html
# either call find_package() or set PYBIND11_FINDPYTHON to ON
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
#pybind11_add_module(${DD2MD} dd2md.cpp) # not working for cmake 3.16
Python_add_library(${DD2MD} dd2md.cpp) # ???same source as above???
target_link_libraries(${DD2MD} PUBLIC pybind11::headers Eigen3::Eigen pybind11::module ${DDOMP_Dir})
set_target_properties(${DD2MD} PROPERTIES
                     CXX_VISIBILITY_PRESET default
                     VISIBILITY_INLINES_HIDDEN ON
                     LIBRARY_OUTPUT_NAME ${DD2MD_OUTPUT_NAME}
                     PREFIX "")
   #INTERPROCEDURAL_OPTIMIZATION ON

#add_library(${DD2MD} MODULE dd2md.cpp)
set_target_properties(${DD2MD} PROPERTIES LIBRARY_OUTPUT_NAME ${DD2MD_OUTPUT_NAME} PREFIX "")

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

# append to DD2MD_LINK_LIBRARIES_LIST: Eigen3::Eigen pybind11::module
list(APPEND DD2MD_LINK_LIBRARIES_LIST Eigen3::Eigen pybind11::module)
if(USE_NOISE_GENERATOR)
   list(APPEND DD2MD_LINK_LIBRARIES_LIST ${FFTW3_LIBRARIES})
endif()
message(STATUS "DD2MD_LINK_LIBRARIES_LIST ${DD2MD_LINK_LIBRARIES_LIST}")

target_link_libraries(${DD2MD} PUBLIC ${DD2MD_LINK_LIBRARIES_LIST})

target_include_directories(${DD2MD} PUBLIC ${DD2MD_INCLUDE_Dir})
target_include_directories(${DD2MD} PUBLIC Eigen3::Eigen)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/Utilities/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/IO/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/Math/)
#target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/MPI/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/Geometry/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/Mesh/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/Lattices/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/FEM/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/GlidePlanes/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/LoopNetwork/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/Quadrature/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/ParticleInteraction/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DiscreteCrackMechanics/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/PolycrystallineMaterials/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DislocationDynamicsIO/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DislocationMicrostructure/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DislocationMobilities/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DislocationDynamicsBase/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DislocationDynamics/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DislocationDynamics/Polycrystals/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DislocationDynamics/Visualization/)
target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/DislocationDynamics/ElasticFields/)
#target_include_directories(${DD2MD} PUBLIC ${MoDELibRootDir}/include/vtk/)

add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
add_subdirectory(${MoDELibRootDir}/src/GlidePlanes GlidePlanes)
add_subdirectory(${MoDELibRootDir}/src/DislocationMobilities DDDmobilities)
add_subdirectory(${MoDELibRootDir}/src/DislocationDynamicsIO DDDIO)
add_subdirectory(${MoDELibRootDir}/src/DislocationMicrostructure DDDmicrostructure)
add_subdirectory(${MoDELibRootDir}/src/PolycrystallineMaterials PolyMaterials)
add_subdirectory(${MoDELibRootDir}/src/DislocationDynamics DDD)
#add_subdirectory(${MoDELibRootDir}/src/vtk VTK)

