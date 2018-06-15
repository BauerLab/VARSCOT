cmake_minimum_required (VERSION 3.0.0)

project (VARSCOT_variants CXX)

set (CMAKE_BUILD_TYPE "Release")

set (CMAKE_PREFIX_PATH "../lib/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "../lib/seqan/include")

find_package (OpenMP REQUIRED)
find_package (SeqAn REQUIRED)
find_package (ZLIB REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

set (CMAKE_CXX_FLAGS "-std=c++14 -fopenmp -O3 ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

add_executable (fasta_writer fasta_writer.cpp)
target_link_libraries (fasta_writer ${SEQAN_LIBRARIES})

add_executable (vcf_loader vcf_loader.cpp)
target_link_libraries (vcf_loader ${SEQAN_LIBRARIES})

add_executable (bam_merger bam_merger.cpp)
target_link_libraries (bam_merger ${SEQAN_LIBRARIES})

add_executable (bam_merger_ref_only bam_merger_ref_only.cpp)
target_link_libraries (bam_merger_ref_only ${SEQAN_LIBRARIES})
