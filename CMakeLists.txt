cmake_minimum_required(VERSION 3.1.0)

project(Metaphor C CXX)

set(META_VERSION_MAJOR 1)
set(META_VERSION_MINOR 0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(META_LIB_NAME metaphor${META_VERSION_MAJOR}.${META_VERSION_MINOR})


file(GLOB META_CORE_SRC
"Core/*.h"
"Core/Include/*.h"
"Core/Src/*.cpp"
)

add_library(
    ${META_LIB_NAME} STATIC
    ${META_CORE_SRC}
)

set_property(TARGET ${META_LIB_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${META_LIB_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

target_include_directories(
    ${META_LIB_NAME} PUBLIC
    Core
)
