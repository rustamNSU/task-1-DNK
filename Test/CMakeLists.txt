cmake_minimum_required(VERSION 3.17)

project (Test)

file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Compile executable from source files in SOURCES variable
add_executable(DNKTest ${SOURCES})
target_link_libraries(DNKTest googletest)


#install (TARGETS LayeredElasticityMatrixUnitTest DESTINATION bin)
