include_directories(${MoDELibRootDir}/include/DislocationDynamicsBase)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials)
include_directories(${MoDELibRootDir}/include/Mesh)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)

target_sources(${PROJECT_NAME} PRIVATE 
   ${CMAKE_CURRENT_SOURCE_DIR}/DislocationDynamicsBase.cpp
   )
