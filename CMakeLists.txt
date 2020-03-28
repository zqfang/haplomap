cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(haplomap)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
#find_package(GSL REQUIRED)

# For profiler on Ubuntu
if (CMAKE_BUILD_TYPE STREQUAL Debug AND NOT MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -ggdb -O2 -m64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ggdb -O2 -m64")
endif ()

# GCC need to link against stdc++fs
if(MSVC)
    set(STDCXXFS "")
else()
    set(STDCXXFS "stdc++fs")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ggdb -O2 -m64")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# add GSL
set(GSL_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/gsl/include)
set(GSL_LIBS ${CMAKE_CURRENT_SOURCE_DIR}/gsl/lib)

# add srcs, headers 
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${GSL_INCLUDE})
link_directories(${GSL_LIBS}) # so, libgsl.so libgslcblas.so is found in ${LIBS}
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src HALPOMAP_CXX) # so, all .cpp file could be found by targets

## targets
# eblocks
add_executable(eblocks ${HALPOMAP_CXX} ${CMAKE_CURRENT_SOURCE_DIR}/ehaploblocks.cpp)
target_include_directories(eblocks PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE})
#target_link_libraries(eblocks )

# ghmap
add_executable(ghmap ${HALPOMAP_CXX} ${CMAKE_CURRENT_SOURCE_DIR}/quantTraitMap.cpp)
target_include_directories(ghmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE})
target_link_libraries(ghmap PUBLIC gsl gslcblas ) # link libgsl.so libgslcblas.so 




