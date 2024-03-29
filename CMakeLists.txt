# 3.16: Oldest CMake version present on Ubuntu LTS
cmake_minimum_required(VERSION 3.16)

# set the project name
project(ESSColor)
project(ESSColor VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# declare libraries
add_subdirectory(lib)
add_subdirectory(src)

    
# add the executable
add_executable(essColorAuxKmerlist src/kmerlist.cpp)
add_executable(essColorAuxMatrixCompress src/compress_matrix.cpp)
add_executable(essColorAuxMatrixDecompress src/decompress_matrix.cpp)
add_executable(genmatrix src/genMatrix.cpp)
add_executable(essColorAuxConfigWriter src/config_writer.cpp)

add_executable(essColorAuxKmersDecompress src/essCompress/decoder.cpp)
add_executable(essColorAuxKmersCompress src/essCompress/ess.cpp)



 



# link libraries
target_link_libraries(essColorAuxMatrixCompress cmph_lib sdsl)
target_link_libraries(essColorAuxMatrixDecompress cmph_lib sdsl)
target_link_libraries(genmatrix essc_lib cxxopts sshash_static)
#target_link_libraries(genmatrix cxxopts sshash_static)