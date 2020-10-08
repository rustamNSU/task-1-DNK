cmake_minimum_required(VERSION 3.17)
project(DNK)

set(CMAKE_CXX_STANDARD 20)

include(ExternalProject)
include(${CMAKE_CURRENT_SOURCE_DIR}/externalGoogleTest.cmake)
add_subdirectory(Test)

#file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
#file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
##message(STATUS ${SOURCES})
#
#add_executable(DNK ${SOURCES})
#target_link_libraries(DNK googletest)