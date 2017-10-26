cmake_minimum_required (VERSION 3.0.0)

project (VARSCOT_search CXX)

set (CMAKE_BUILD_TYPE "Release")
set (CXX11_FOUND "1")

set (CMAKE_PREFIX_PATH "../lib/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "../lib/seqan/include")

find_package (OpenMP REQUIRED)
find_package (SeqAn REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

# -DENABLE_DEBUG_MACRO builds alignments during search to help debugging search schemes
# This is only activated in Debug mode since it increases the running time.
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DENABLE_DEBUG_MACRO -DSEQAN_DISABLE_VERSION_CHECK=YES")
set (CMAKE_CXX_FLAGS "-std=c++14 -fopenmp ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -DSEQAN_DISABLE_VERSION_CHECK=YES")

add_executable (bidir_index bidir_index.cpp)
target_link_libraries (bidir_index ${SEQAN_LIBRARIES})

add_executable (bidir_mapping bidir_mapping.cpp)
target_link_libraries (bidir_mapping ${SEQAN_LIBRARIES})
