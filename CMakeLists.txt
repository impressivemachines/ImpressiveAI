cmake_minimum_required(VERSION 2.8)

project(Metaphor C CXX)

set (META_VERSION_MAJOR 0)
set (META_VERSION_MINOR 1)

file(GLOB META_SRC_CORE     Src/Core/*.cpp)
file(GLOB META_SRC_IO       Src/IO/*.cpp)
file(GLOB META_SRC_Sound    Src/Sound/*.cpp)
file(GLOB META_SRC_Vision   Src/Vision/*.cpp)

add_library(metaphor STATIC test.c)