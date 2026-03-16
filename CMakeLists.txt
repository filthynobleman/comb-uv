cmake_minimum_required(VERSION 3.23.0)
project(Diskify LANGUAGES C CXX)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/lib")
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/lib")
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/bin")
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

# Dependencies
include_directories("${CMAKE_SOURCE_DIR}/ext/eigen")
include_directories("${CMAKE_SOURCE_DIR}/ext/libigl/include")
include_directories("${CMAKE_SOURCE_DIR}/ext/stb")

add_library(STB STATIC "${CMAKE_SOURCE_DIR}/src/stb/stb.cpp")

# Project headers
include_directories("${CMAKE_SOURCE_DIR}/include")

# Parallelism
find_package(OpenMP)


# Project binaries
add_library(DFY STATIC  "${CMAKE_SOURCE_DIR}/src/dfy/mesh.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/manifoldmesh.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/graph.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/sampler.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/region.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/segmentation.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/embedding.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/tutte.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/harmonic.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/arap.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/conformal.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/minheap.cpp"
                        "${CMAKE_SOURCE_DIR}/src/dfy/io.cpp")
target_link_libraries(DFY STB)
if (OpenMP_FOUND)
    target_link_libraries(DFY OpenMP::OpenMP_CXX)
endif()
set_target_properties(DFY PROPERTIES CXX_STANDARD 17)


# Applications
add_executable(UVPatches "${CMAKE_SOURCE_DIR}/src/apps/uvpatches.cpp")
target_link_libraries(UVPatches DFY)
set_target_properties(UVPatches PROPERTIES CXX_STANDARD 17)


add_executable(Diskify "${CMAKE_SOURCE_DIR}/src/apps/diskify.cpp")
target_link_libraries(Diskify DFY)
set_target_properties(Diskify PROPERTIES CXX_STANDARD 17)


add_executable(GraphicalAbstract "${CMAKE_SOURCE_DIR}/src/apps/graphicalabstract.cpp")
target_link_libraries(GraphicalAbstract DFY)
set_target_properties(GraphicalAbstract PROPERTIES CXX_STANDARD 17)


add_executable(VoronoiRefinementFigure "${CMAKE_SOURCE_DIR}/src/apps/vororefine.cpp")
target_link_libraries(VoronoiRefinementFigure DFY)
set_target_properties(VoronoiRefinementFigure PROPERTIES CXX_STANDARD 17)