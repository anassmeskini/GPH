cmake_minimum_required (VERSION 3.0.0)
project(heuristic VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}-Wall")

file(GLOB SOURCES "src/*.cpp")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_executable(heuristic ${SOURCES})

find_package(Boost COMPONENTS iostreams REQUIRED)
find_package(TBB REQUIRED tbb)
find_library(Cplex REQUIRED)

add_compile_definitions(CONCERT_CPLEX_FOUND)

add_compile_definitions(IL_STD)

include_directories(${PROJECT_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS} ${CPLEX_CONCERT_INCLUDE_DIR} ${CPLEX_INCLUDE_DIR})

target_link_libraries(heuristic PRIVATE ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${CPLEX_CONCERT_LIBRARY} ${CPLEX_ILOCPLEX_LIBRARY} dl ${CPLEX_LIBRARY} pthread)
