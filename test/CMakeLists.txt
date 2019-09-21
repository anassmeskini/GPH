
add_executable(unit_test "${CMAKE_SOURCE_DIR}/test/test.cpp" ${SOURCES})

target_link_libraries(unit_test PRIVATE 
                      ${TBB_IMPORTED_TARGETS}
                      ${LPSOLVER_LIBRARIES})

add_test(unit_test test)
