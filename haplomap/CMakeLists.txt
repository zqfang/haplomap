project(haplomap)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

# For profiler on Ubuntu
if (CMAKE_BUILD_TYPE STREQUAL Debug AND NOT MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -g -O2 -m64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O2 -m64")
endif ()

# GCC need to link against stdc++fs
#if(MSVC)
#    set(STDCXXFS "")
#else()
#    set(STDCXXFS "stdc++fs")
#endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# add srcs, headers 
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src HAPLOMAP_SRC)

add_executable(haplomap main.cpp "${HAPLOMAP_SRC}")
target_include_directories(haplomap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS})
target_link_libraries(haplomap PUBLIC ${GSL_LIBRARIES}) # link libgsl.so libgslcblas.so
#install binary to directory you specified
install(TARGETS haplomap 
        CONFIGURATIONS Release
        RUNTIME DESTINATION "${HAPLOMAP_INSTALL_BINDIR}")
      #endif()







