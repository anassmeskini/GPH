add_executable(gph "main.cpp")

target_include_directories(gph PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/external)
target_include_directories(gph SYSTEM PRIVATE ${LPSOLVER_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
target_link_libraries(gph PRIVATE mylib)
