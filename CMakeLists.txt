cmake_minimum_required(VERSION 3.5)
project(HBCGM)
set(CMAKE_CXX_STANDARD 11)

# # Download and unpack googletest at configure time
# configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
# execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
# execute_process(COMMAND ${CMAKE_COMMAND} --build .
#         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )

# # Prevent GoogleTest from overriding our compiler/linker options
# # when building with Visual Studio
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# # Add googletest directly to our build. This adds
# # the following targets: gtest, gtest_main, gmock
# # and gmock_main
# add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
#         ${CMAKE_BINARY_DIR}/googletest-build)

# # The gtest/gmock targets carry header search path
# # dependencies automatically when using CMake 2.8.11 or
# # later. Otherwise we have to add them here ourselves.
# if (CMAKE_VERSION VERSION_LESS 2.8.11)
#     include_directories("${gtest_SOURCE_DIR}/include" "${gmock_SOURCE_DIR}/include")
# endif()
# # Now, you could use unitest inside u projects
# enable_testing()

#find_package(GSL REQUIRED)
# define some globals
SET(SOURCECODE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# output dir
# add GSL
set(GSL_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/gsl/include)
set(GSL_LIBS ${CMAKE_CURRENT_SOURCE_DIR}/gsl/lib)
include_directories(${GSL_INCLUDE})
link_directories(${GSL_LIBS})  # so, libgsl.so libgslcblas.so are found in ${GSL_LIBS}
set(GSL_LIBRARIES gsl gslcblas)

file(GLOB HBCGM_SRC ${CMAKE_SOURCE_DIR}/*/src/*.cpp)

# Set the install prefix
if(DEFINED CMAKE_INSTALL_PREFIX)
   set(HAPLOMAP_INSTALL_BINDIR "${CMAKE_INSTALL_PREFIX}" )
else()
   set(HAPLOMAP_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}")
endif()


########## set BUID Flag #########
option(ENABLE_HAPLOMAP_BUILD
       "Build the haplomap"
       ON
)

option(ENABLE_UNITEST_BUILD 
       "Build the UNITEST"
       OFF
)

if(ENABLE_HAPLOMAP_BUILD)
    add_subdirectory(haplomap)
else(ENABLE_HAPLOMAP_BUILD)
    message(STATUS "Skipping haplomap build")
endif(ENABLE_HAPLOMAP_BUILD)

if(ENABLE_UNITEST_BUILD)
    add_subdirectory(test)
else(ENABLE_UNITEST_BUILD)
    message(STATUS "Skipping unitest build")
endif(ENABLE_UNITEST_BUILD)

