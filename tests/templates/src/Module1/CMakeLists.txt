message(STATUS "adding sources Module1")
message(STATUS "including ${MoDELibRootDir}/include/Module1/")
include_directories(${MoDELibRootDir}/include/Module1/)
add_library(Module1 STATIC MyClassTemplate.cpp)

