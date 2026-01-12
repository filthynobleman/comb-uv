cmake_minimum_required(VERSION 3.23.0)
project(Diskify LANGUAGES C CXX)


# Dependencies
include_directories("${CMAKE_SOURCE_DIR}/ext/eigen")
include_directories("${CMAKE_SOURCE_DIR}/ext/libigl/include")

# Project headers
include_directories("${CMAKE_SOURCE_DIR}/include")


# Project binaries
add_library(DFY STATIC  "${CMAKE_SOURCE_DIR}/src/dfy/mesh.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/manifoldmesh.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/graph.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/sampler.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/region.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/segmentation.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/minheap.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/io.cpp")
set_target_properties(DFY PROPERTIES CXX_STANDARD 17)


# Application
add_executable(Test "${CMAKE_SOURCE_DIR}/src/main.cpp")
target_link_libraries(Test DFY)
set_target_properties(Test PROPERTIES CXX_STANDARD 17)