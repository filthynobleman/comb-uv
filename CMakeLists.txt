cmake_minimum_required(VERSION 3.23.0)
project(Diskify LANGUAGES C CXX)


# Dependencies
include_directories("${CMAKE_SOURCE_DIR}/ext/eigen")
include_directories("${CMAKE_SOURCE_DIR}/ext/libigl/include")

# Project headers
include_directories("${CMAKE_SOURCE_DIR}/include")