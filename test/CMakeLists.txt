project(test)


set(Pleiades_INC ${CMAKE_SOURCE_DIR}/haplomap/include)

## math files
# file(GLOB Pleiades_SRC ${CMAKE_SOURCE_DIR}/haplomap/src/*.cpp)
include_directories(${Pleiades_INC})
# Find all source files in a directory
aux_source_directory(${CMAKE_SOURCE_DIR}/haplomap/src SRC) # so, all .cpp file could be found by targets
aux_source_directory(${PROJECT_SOURCE_DIR} TEST_SRC)
# First Simple Google Test
add_executable(unitest unitest.cpp ${SRC} ghmap.test.cpp)# ${TEST_SRC})
#target_include_directories(unitest ${Pleiades_INC})
target_link_libraries(unitest gtest gmock_main ${GSL_LIBRARIES})
add_test(NAME Pleiades_test COMMAND unitest)
