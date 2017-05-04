include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/common")

add_executable(point_test point_test.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:distance> $<TARGET_OBJECTS:matrix>)
target_link_libraries(point_test gtest_main)
add_test(NAME point_test COMMAND point_test)
