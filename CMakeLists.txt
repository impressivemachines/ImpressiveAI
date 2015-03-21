cmake_minimum_required(VERSION 2.8.11)

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

#file(GLOB META_VISION_SRC
#"Vision/*.h"
#"Vision/Include/*.h"
#"Vision/Src/*.cpp"
#)

add_library(
    ${META_LIB_NAME} STATIC
    ${META_CORE_SRC}
#    ${META_VISION_SRC}
)

target_include_directories(
    ${META_LIB_NAME} PUBLIC
    Core
#    Vision
)

#add_executable(
#    debughacks
#    Test/Src/debughacks.cpp
#)

#target_link_libraries(
#    debughacks
#    ${META_LIB_NAME}
#)
