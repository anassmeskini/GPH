cmake_minimum_required (VERSION 3.0.0)
project(GPH LANGUAGES CXX
            DESCRIPTION
            "General-purpose heuristics for mixed-integer programming")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB SOURCES "src/core/*.cpp" 
                  "src/interfaces/*.cpp"  
                  "src/io/*.cpp" 
                  "src/methods/*.cpp")

add_executable(gph "src/main.cpp" ${SOURCES})

add_executable(scp "example/SetCover.cpp" "example/GreedyHeuristic.cpp" ${SOURCES})

#compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# dependencies
find_package(TBB REQUIRED tbb tbbmalloc_proxy)

# use the header only version of fmtlib
add_compile_definitions(FMT_HEADER_ONLY)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# solvers
if(SOLVER STREQUAL "cplex")
   find_package(Cplex REQUIRED)
elseif(SOLVER STREQUAL "soplex")
   find_package(SOPLEX CONFIG)
elseif(SOLVER STREQUAL "glpk")
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
   message(FATAL_ERROR "LP solver not found")
endif()

target_include_directories(gph PRIVATE ${PROJECT_SOURCE_DIR}/src
                    ${PROJECT_SOURCE_DIR}/external
                    SYSTEM
                    ${TBB_INCLUDE_DIRS}
                    ${LPSOLVER_INCLUDE_DIRS})

target_include_directories(scp PRIVATE ${PROJECT_SOURCE_DIR}/src
                    ${PROJECT_SOURCE_DIR}/external
                    SYSTEM
                    ${TBB_INCLUDE_DIRS}
                    ${LPSOLVER_INCLUDE_DIRS})

target_link_libraries(gph PRIVATE 
                      TBB::tbb
                      TBB::tbbmalloc_proxy
                      ${LPSOLVER_LIBRARIES} m stdc++fs)

target_link_libraries(scp PRIVATE
                      TBB::tbb
                      TBB::tbbmalloc_proxy
                      ${LPSOLVER_LIBRARIES} m stdc++fs)

add_subdirectory(test)
