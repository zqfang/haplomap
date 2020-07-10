project(test)

## math files
file(GLOB HBCGM_SOURCES "*.cpp" "*.c" "${CMAKE_SOURCE_DIR}/*/src/*.cpp" )
file(GLOB HBCGM_INCLUDES "${CMAKE_SOURCE_DIR}/*/include/" )

include_directories(${HBCGM_INCLUDES})
# Find all source files in a directory
# aux_source_directory(${CMAKE_SOURCE_DIR}/haplomap/src SRC) # so, all .cpp file could be found by targets
# aux_source_directory(${PROJECT_SOURCE_DIR} TEST_SRC)
# First Simple Google Test
add_executable(unitest ${HBCGM_SOURCES})# ${TEST_SRC})
#target_include_directories(unitest ${Pleiades_INC})
target_link_libraries(unitest gtest gmock_main ${GSL_LIBRARIES})
add_test(NAME HBCGM_test COMMAND unitest)
