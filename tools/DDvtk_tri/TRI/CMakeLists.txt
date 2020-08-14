cmake_minimum_required(VERSION 3.0)

project(tri CXX)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")

include_directories(include)

add_definitions(-DTRILIBRARY)
add_definitions(-DANSI_DECLARATORS)
add_definitions(-DCDT_ONLY)



file(GLOB src_files "src/assert.cpp"
					"src/del_impl.cpp")


add_executable(tri src/main.cpp ${src_files})
