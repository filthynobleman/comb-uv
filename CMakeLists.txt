cmake_minimum_required(VERSION 3.23.0)
project(Diskify LANGUAGES C CXX)


# Dependencies
include_directories("${CMAKE_SOURCE_DIR}/ext/eigen")
include_directories("${CMAKE_SOURCE_DIR}/ext/libigl/include")
include_directories("${CMAKE_SOURCE_DIR}/ext/stb")

add_library(STB STATIC "${CMAKE_SOURCE_DIR}/src/stb/stb.cpp")

# Project headers
include_directories("${CMAKE_SOURCE_DIR}/include")


# Project binaries
add_library(DFY STATIC  "${CMAKE_SOURCE_DIR}/src/dfy/mesh.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/manifoldmesh.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/quadmesh.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/graph.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/sampler.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/region.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/segmentation.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/embedding.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/tutte.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/harmonic.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/arap.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/conformal.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/gimage.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/imgsampler.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/minheap.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/io.cpp")
target_link_libraries(DFY STB)
set_target_properties(DFY PROPERTIES CXX_STANDARD 17)


# Application
add_executable(Test "${CMAKE_SOURCE_DIR}/src/main.cpp")
target_link_libraries(Test DFY)
set_target_properties(Test PROPERTIES CXX_STANDARD 17)


add_executable(UVPatches "${CMAKE_SOURCE_DIR}/src/apps/uvpatches.cpp")
target_link_libraries(UVPatches DFY)
set_target_properties(UVPatches PROPERTIES CXX_STANDARD 17)