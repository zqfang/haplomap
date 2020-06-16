  
include_directories("${PROJECT_SOURCE_DIR}")

include_directories(${CMAKE_SOURCE_DIR}/haplomap/include)
include_directories(${GSL_INCLUDE})
link_directories(${GSL_LIBS}) # so, libgsl.so libgslcblas.so is found in ${LIBS}
aux_source_directory(${CMAKE_SOURCE_DIR}/haplomap/src HAPLOMAP_SRC) # so, all .cpp file could be found by targets





add_executable(main main.cpp)
target_link_libraries(main gtest)
# First Simple Google Test
add_executable(simple unitest.cpp ${HAPLOMAP_SRC})
target_link_libraries(simple gtest gmock_main)
add_test(NAME simple_test COMMAND simple)