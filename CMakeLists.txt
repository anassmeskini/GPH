cmake_minimum_required (VERSION 3.0.0)
project(LPS VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

file(GLOB SOURCES "src/main.cpp" "src/core/*.cpp" "src/interface/*.cpp")
add_executable(lps ${SOURCES})

find_package(Boost COMPONENTS iostreams REQUIRED)
find_package(TBB REQUIRED tbb)

#use the header only version of fmtlib
add_compile_definitions(FMT_HEADER_ONLY)

#solvers
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(SOLVER STREQUAL "cpx")
   find_package(Cplex REQUIRED)
elseif(SOLVER STREQUAL "spx")
   find_package(SOPLEX CONFIG)
elseif(SOLVER STREQUAL "glp")
   find_package(GLPK)
elseif(SOLVER STREQUAL "")
   message(FATAL_ERROR "LP solver not specified")
else()
   message(FATAL_ERROR "LP solver unknown")
endif()

if(CPLEX_FOUND)
   set(LPSOLVER_INCLUDE_DIRS ${CPLEX_CONCERT_INCLUDE_DIR} ${CPLEX_INCLUDE_DIR})
   set(LPSOLVER_LIBRARIES ${CPLEX_CONCERT_LIBRARY} ${CPLEX_ILOCPLEX_LIBRARY} dl ${CPLEX_LIBRARY} pthread)
   add_compile_definitions(CONCERT_CPLEX_FOUND)
   add_compile_definitions(IL_STD)
elseif(SOPLEX_FOUND)
   set(LPSOLVER_INCLUDE_DIRS ${SOPLEX_INCLUDE_DIRS})
   set(LPSOLVER_LIBRARIES ${SOPLEX_LIBRARIES})
   add_compile_definitions(SOPLEX_FOUND)
elseif(GLPK_FOUND)
   set(LPSOLVER_INCLUDE_DIRS ${GLPK_INCLUDE_DIR})
   set(LPSOLVER_LIBRARIES ${GLPK_LIBRARY})
   add_compile_definitions(GLPK_FOUND)
else()
   message(FATAL_ERROR "LP solver NOT found")
endif()

include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/external SYSTEM ${Boost_INCLUDE_DIRS} ${LPSOLVER_INCLUDE_DIRS})
target_link_libraries(lps PRIVATE ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${LPSOLVER_LIBRARIES})
