cmake_minimum_required (VERSION 3.0.0)

project(LPS VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

file(GLOB SOURCES "src/*.cpp")
add_executable(lps ${SOURCES})

find_package(Boost COMPONENTS iostreams REQUIRED)
find_package(TBB REQUIRED tbb)

#solvers
if(LPS STREQUAL "cpx")
   find_package(Cplex)
   if(CPLEX_FOUND)
      set(LPSOLVER_INCLUDE_DIRS ${CPLEX_CONCERT_INCLUDE_DIR} ${CPLEX_INCLUDE_DIR})
      set(LPSOLVER_LIBRARIES ${CPLEX_CONCERT_LIBRARY} ${CPLEX_ILOCPLEX_LIBRARY} dl ${CPLEX_LIBRARY} pthread)
      add_compile_definitions(CONCERT_CPLEX_FOUND)
      add_compile_definitions(IL_STD)
   else()
      message(FATAL_ERROR Cplex not found)
   endif()
elseif(LPS STREQUAL "spx")
   find_package(SOPLEX CONFIG)
   if(SOPLEX_FOUND)
      add_compile_definitions(SOPLEX_FOUND)
      set(LPSOLVER_INCLUDE_DIRS ${SOPLEX_INCLUDE_DIRS})
      set(LPSOLVER_LIBRARIES ${SOPLEX_LIBRARIES})
   else()
      message(FATAL_ERROR Soplex not found)
   endif()
else()
   message(FATAL_ERROR LP solver not specified)
endif()

# SYSTEM disables warnings from external headers
include_directories(${PROJECT_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS} SYSTEM ${LPSOLVER_INCLUDE_DIRS})

target_link_libraries(lps PRIVATE ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${LPSOLVER_LIBRARIES})

target_compile_options(lps PRIVATE "-Wall")
