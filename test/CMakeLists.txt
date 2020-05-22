  
include_directories("${PROJECT_SOURCE_DIR}")

add_executable(main main.cpp)
target_link_libraries(main gtest)

add_executable(example example.cpp unitest.cpp)
target_link_libraries(example gtest gmock_main)
add_test(NAME example_test COMMAND example)


# First Simple Google Test
add_executable(simple unitest.cpp)
target_link_libraries(simple gtest gmock_main)
add_test(NAME simple_test COMMAND simple)