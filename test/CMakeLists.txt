include_directories("${PROJECT_SOURCE_DIR}")

# First Simple Google Test
add_executable(simple unitest.cpp)
target_link_libraries(simple gtest gmock_main)
add_test(NAME simple_test COMMAND simple)