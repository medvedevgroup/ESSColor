# 3.16: Oldest CMake version present on Ubuntu LTS
cmake_minimum_required(VERSION 3.16)

# set the project name
project(ESSColor)
project(ESSColor VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# declare libraries
add_subdirectory(lib)
add_subdirectory(src)

    
# add the executable
add_executable(kmerlist src/kmerlist.cpp)
add_executable(compress src/compress.cpp)
add_executable(decompress src/decompress.cpp)
add_executable(genmatrix src/genMatrix.cpp)
add_executable(amatur src/amatur.cpp)


# link libraries
target_link_libraries(compress cmph_lib sdsl)
target_link_libraries(decompress cmph_lib sdsl)
target_link_libraries(genmatrix essc_lib cxxopts)
target_link_libraries(amatur sshash)
